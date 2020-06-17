import re
import json
import uuid
import boto3
from datetime import datetime, timezone, timedelta
from botocore.client import Config


S3_BUCKET = 'codfreq-assets.hivdb.org'
S3 = boto3.client('s3', config=Config(signature_version='s3v4'))
EXPIRATION = timedelta(days=7)


def utcnow_text():
    return datetime.now(tz=timezone.utc).isoformat()


def utcnow():
    return datetime.now(tz=timezone.utc)


def parse_dt(value):
    return datetime.strptime(value, '%Y-%m-%dT%H:%M:%S.%f%z')


def load_taskmeta(uniqkey):
    obj = S3.get_object(
        Bucket=S3_BUCKET,
        Key='tasks/{}/taskmeta.json'.format(uniqkey)
    )
    return json.load(obj['Body'])


def is_expired(data):
    last_updated = data['lastUpdatedAt'] = parse_dt(data['lastUpdatedAt'])
    return utcnow() - EXPIRATION > last_updated


def can_upload(data):
    return data['status'] == 'created'


def can_trigger(data):
    return data['status'] == 'wait-client-upload'


def save_taskmeta(uniqkey, status, payload=None):
    now = utcnow_text()
    S3.put_object(
        Bucket=S3_BUCKET,
        Key='tasks/{}/taskmeta.json'.format(uniqkey),
        ContentType='application/json',
        Body=json.dumps({
            'taskKey': uniqkey,
            'lastUpdatedAt': now,
            'status': status,
            'payload': payload
        }).encode('U8')
    )
    return now


def create_task(request):
    """Step 1: create task meta file"""
    uniqkey = str(uuid.uuid4())
    now = save_taskmeta(uniqkey, 'created')
    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': now,
        'status': 'created'
    }, 200


def trigger_runner(request):
    """Step 3: trigger codfreq-runner"""
    uniqkey = request.get('taskKey')
    if not uniqkey:
        return {"error": "'taskKey' not provided"}, 400
    try:
        taskmeta = load_taskmeta(uniqkey)
    except S3.exceptions.NoSuchKey:
        return {"error": "taskKey not found"}, 404
    if is_expired(taskmeta):
        return {"error": "this task is expired"}, 404
    if not can_trigger(taskmeta):
        return {"error": "this task is not ready for a runner"}, 400
    payload = taskmeta['payload']
    fastq_files = set(payload['fastqFiles'])
    result = S3.list_objects_v2(
        Bucket=S3_BUCKET,
        MaxKeys=len(fastq_files) + 1,
        Prefix='tasks/{}/'.format(uniqkey)
    )
    if result['IsTruncated']:
        return {"error": "too many files are found for this task"}, 400
    known_files = {f['Key'].rsplit('/', 1)[-1] for f in result['Contents']}
    known_files.remove('taskmeta.json')
    missing_files = fastq_files - known_files
    if missing_files:
        return {'error': 'not all fastq files are uploaded',
                'missingFiles': sorted(missing_files)}, 400
    return {'ok': "let's go!"}, 200


def direct_upload(request):
    """Step 2: retrieve FASTQ files uploading credential"""
    uniqkey = request.get('taskKey')
    filenames = request.get('fileNames')
    if not uniqkey:
        return {"error": "'taskKey' not provided"}, 400

    if not filenames:
        return {"error": "'fileNames' not provided"}, 400

    if not isinstance(filenames, list):
        return {"error": "'fileNames' must be a list of strings"}, 400

    for filename in filenames:
        fnlower = filename.lower()
        if not fnlower.endswith('.fastq') and \
                not fnlower.endswith('.fastq.gz'):
            return {
                'error': 'attempt to upload a non-FASTQ file {!r}'
                .format(filename)
            }, 400
        matches = re.findall(r'[^\w.()-]', filename)
        matches = ''.join(sorted(set(matches)))
        if matches:
            return {
                'error': 'filename {!r} contains invalid character(s): {!r}'
                .format(filename, matches)
            }, 400

    try:
        taskmeta = load_taskmeta(uniqkey)
    except S3.exceptions.NoSuchKey:
        return {"error": "taskKey not found"}, 404
    if is_expired(taskmeta):
        return {"error": "this task is expired"}, 404
    if not can_upload(taskmeta):
        return {"error": "this task is no longer accept new files"}, 404

    result = {
        'taskKey': uniqkey,
        'status': 'wait-client-upload'
    }
    sigs = []
    for filename in filenames:
        sigs.append(
            S3.generate_presigned_post(
                Bucket=S3_BUCKET,
                Key='tasks/{}/{}'.format(uniqkey, filename),
                Fields={'content-type': 'binary/octet-stream'},
                Conditions=[{'content-type': 'binary/octet-stream'}],
                ExpiresIn=3600
            )
        )
    result['presignedPosts'] = sigs
    now = save_taskmeta(uniqkey, 'wait-client-upload', {
        'fastqFiles': filenames
    })
    result['lastUpdatedAt'] = now
    return result, 200


def method_not_found(request):
    return {"error": "method not found"}, 404


def dispatch(event, context):
    handlername = event['pathParameters']['handler']
    request = json.loads(event['body'] or 'null')
    if handlername == 'create-task':
        method = create_task
    elif handlername == 'direct-upload':
        method = direct_upload
    elif handlername == 'trigger-runner':
        method = trigger_runner
    else:
        method = method_not_found
    response_body, status = method(request)

    return {
        'statusCode': status,
        'headers': {
            'Access-Control-Allow-Origin': '*',
            'Content-Type': 'application/json',
        },
        'body': json.dumps(response_body)
    }
