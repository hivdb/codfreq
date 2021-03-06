import re
import csv
import json
import uuid
import boto3
from io import StringIO
from datetime import datetime, timezone, timedelta
from botocore.client import Config


ECS_REGION = 'us-west-2'
ECS_CLUSTER = 'codfreq-runner'
ECS_TASK_DEFINITION = 'codfreq-runner'
ECS_NETWORK_CONFIG = {
    'awsvpcConfiguration': {
        'subnets': [
            'subnet-02f7aeefbe32139e6',
            'subnet-0d99ba6fe930466bc'
        ],
        'securityGroups': [
            'sg-0ab407848018a5604'
        ],
        'assignPublicIp': 'ENABLED'
    }
}
EXPIRATION = timedelta(days=7)


S3_BUCKET = 'codfreq-assets.hivdb.org'
S3 = boto3.client('s3', config=Config(signature_version='s3v4'))
ECS = boto3.client('ecs', region_name=ECS_REGION)
LOGS = boto3.client('logs', region_name=ECS_REGION)


def utcnow_text():
    return datetime.now(tz=timezone.utc).isoformat()


def utcnow():
    return datetime.now(tz=timezone.utc)


def parse_dt(value):
    return datetime.strptime(value, '%Y-%m-%dT%H:%M:%S.%f%z')


def is_expired(data):
    last_updated = data['lastUpdatedAt']
    return utcnow() - EXPIRATION > last_updated


def can_upload(data):
    return data['status'] == 'created'


def can_trigger(data):
    return data['status'] == 'await-client-upload'


def is_triggered(data):
    return data['status'] == 'await-runner'


def is_success(data):
    return data['status'] == 'success'


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


def load_taskmeta(uniqkey):
    obj = S3.get_object(
        Bucket=S3_BUCKET,
        Key='tasks/{}/taskmeta.json'.format(uniqkey)
    )
    data = json.load(obj['Body'])
    data['lastUpdatedAt'] = parse_dt(data['lastUpdatedAt'])
    return data


def verify_refs(refs):
    for ref in refs:
        try:
            S3.head_object(
                Bucket=S3_BUCKET,
                Key='refs/{}'.format(ref)
            )
        except S3.exceptions.NoSuchKey:
            return False
    return True


def create_task(request):
    """Step 1: create task meta file"""
    uniqkey = str(uuid.uuid4())
    now = save_taskmeta(uniqkey, 'created')
    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': now,
        'status': 'created'
    }, 200


def yield_codfreqs(path_prefix):
    kw = dict(
        Bucket=S3_BUCKET,
        Prefix=path_prefix
    )
    while True:
        result = S3.list_objects_v2(**kw)
        for obj in result['Contents']:
            if not obj['Key'].endswith('.codfreq'):
                continue
            yield obj['Key'], S3.get_object(Bucket=S3_BUCKET, Key=obj['Key'])
        next_token = result.get('NextContinuationToken')
        if not next_token:
            break
        kw['ContinuationToken'] = next_token


def fetch_codfreqs(request):
    """Step 5: fetch codfreq files"""
    uniqkey = request.get('taskKey')
    if not uniqkey:
        return {"error": "'taskKey' not provided"}, 400
    try:
        taskmeta = load_taskmeta(uniqkey)
    except S3.exceptions.NoSuchKey:
        return {"error": "taskKey not found"}, 404
    if is_expired(taskmeta):
        return {"error": "this task is expired"}, 404
    if not is_success(taskmeta):
        return {"error": "this task is not finished yet"}, 400
    path_prefix = 'tasks/{}/'.format(uniqkey)
    codfreqs = {}
    for key, obj in yield_codfreqs(path_prefix):
        name = key.rsplit('.', 2)[0]
        name = name.rsplit('/', 1)[-1]
        name = '{}.codfreq'.format(name)
        gpmap = {}
        fp = StringIO(obj['Body'].read().decode('utf-8-sig'))
        for row in csv.DictReader(fp):
            gene = row['gene']
            pos = int(row['position'])
            total = int(row['total'])
            codon = row['codon']
            if len(codon) < 3:
                continue
            count = int(row['count'])
            refaa = row['refaa']
            aa = row['aa']
            percent = float(row['percent'])
            mean_quality_score = float(row['mean_quality_score'])
            if (gene, pos) not in gpmap:
                gpmap[(gene, pos)] = {
                    'gene': gene,
                    'position': pos,
                    'totalReads': total,
                    'refAminoAcid': refaa,
                    'allCodonReads': []
                }
            gpmap[(gene, pos)]['allCodonReads'].append({
                'codon': codon,
                'reads': count,
                'aminoAcid': aa,
                'percent': percent,
                'meanQualityScore': mean_quality_score
            })
        codfreqs.setdefault(name, []).extend(gpmap.values())
    codfreqs = [{
        'name': name,
        'allReads': all_reads
    } for name, all_reads in codfreqs.items()]

    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': taskmeta['lastUpdatedAt'].isoformat(),
        'status': taskmeta['status'],
        'codfreqs': codfreqs
    }, 200


def fetch_runner_logs(request):
    """Step 4: fetch codfreq-runner logs"""
    uniqkey = request.get('taskKey')
    if not uniqkey:
        return {"error": "'taskKey' not provided"}, 400
    try:
        taskmeta = load_taskmeta(uniqkey)
    except S3.exceptions.NoSuchKey:
        return {"error": "taskKey not found"}, 404
    if is_expired(taskmeta):
        return {"error": "this task is expired"}, 404
    if not is_triggered(taskmeta):
        return {"error": "this task is not triggered yet"}, 400
    payload = taskmeta['payload']
    ecs_task_ids = payload['ecsTaskIds']
    num_tasks = len(ecs_task_ids)
    start_time = request.get('startTime')
    if start_time:
        try:
            start_time = [int(n) for n in start_time.split(',')]
        except (TypeError, ValueError):
            return {"error": "Invalid 'startTime'"}, 400
        if len(start_time) != num_tasks:
            return {
                "error": "'startTime' does not match number of ECS tasks"
            }, 400
    else:
        start_time = [1] * num_tasks

    results = []
    limit = 10000
    finished = set(payload.get('finishedEcsTaskIds', []))
    for task_id, stime in zip(ecs_task_ids, start_time):
        events = []
        kw = dict(
            logGroupName='/ecs/{}'.format(ECS_TASK_DEFINITION),
            logStreamName='ecs/{}/{}'.format(ECS_TASK_DEFINITION,
                                             task_id),
            startTime=stime, limit=limit)
        done = False
        while True:
            try:
                data = LOGS.get_log_events(**kw)
            except LOGS.exceptions.ResourceNotFoundException:
                break
            partial = data.get('events', [])
            if not partial:
                break
            for event in partial:
                try:
                    msg = json.loads(event['message'])
                    events.append({
                        **msg,
                        'timestamp': event['timestamp']
                    })
                    if msg['op'] == 'done':
                        done = True
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception:
                    # ignore non-json message
                    pass
            if done or len(partial) < limit:
                break
            kw['nextToken'] = data['nextForwardToken']
        if done:
            finished.add(task_id)
        results.append({
            'ecsTaskId': task_id,
            'events': events,
            'done': done
        })
    status = taskmeta['status']
    if len(finished) == num_tasks:
        status = 'success'
    payload['finishedEcsTaskIds'] = sorted(finished)
    now = save_taskmeta(uniqkey, status, payload)

    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': now,
        'status': status,
        'taskEvents': results
    }, 200


def trigger_runner(request):
    """Step 3: trigger codfreq-runner"""
    uniqkey = request.get('taskKey')
    runners = request.get('runners')
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
    path_prefix = 'tasks/{}/'.format(uniqkey)
    result = S3.list_objects_v2(
        Bucket=S3_BUCKET,
        MaxKeys=len(fastq_files) + 1,
        Prefix=path_prefix
    )
    if result['IsTruncated']:
        return {"error": "too many files are found for this task"}, 400
    known_files = {f['Key'].rsplit('/', 1)[-1] for f in result['Contents']}
    known_files.remove('taskmeta.json')
    missing_files = fastq_files - known_files
    if missing_files:
        return {'error': 'not all fastq files are uploaded',
                'missingFiles': sorted(missing_files)}, 400
    if not runners:
        return {'error': "'runners' is empty"}, 400
    if not isinstance(runners, list):
        return {'error': "invalid format of 'runners'"}, 400

    all_refs = set()
    for runner in runners:
        if not isinstance(runner, dict):
            return {'error': "invalid format of 'runners'"}, 400
        ref = runner.get('reference')
        all_refs.add(ref)
        gene = runner.get('gene')
        if not gene:
            return {
                'error': "invalid format of 'runners': 'gene' not provided"
            }, 400

    if not verify_refs(all_refs):
        return {
            'error': "invalid format of 'runners': unsupport 'reference'"
        }

    task_ids = []
    for runner in runners:
        ref = runner['reference']
        gene = runner['gene']
        result = ECS.run_task(
            platformVersion='1.4.0',
            cluster=ECS_CLUSTER,
            launchType='FARGATE',
            count=1,
            taskDefinition=ECS_TASK_DEFINITION,
            networkConfiguration=ECS_NETWORK_CONFIG,
            overrides={
                'containerOverrides': [{
                    'name': 'codfreq-runner',
                    'command': [
                        'bin/align-all-docker',
                        '-r', ref,
                        '-g', gene,
                        '-p', path_prefix
                    ]
                }]
            }
        )
        task_id = result['tasks'][0]['taskArn'].rsplit(':task/', 1)[-1]
        task_ids.append(task_id)
    payload['ecsTaskIds'] = task_ids
    now = save_taskmeta(uniqkey, 'await-runner', payload)
    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': now,
        'status': 'await-runner'
    }, 200


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
        'status': 'await-client-upload'
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
    now = save_taskmeta(uniqkey, 'await-client-upload', {
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
    elif handlername == 'fetch-runner-logs':
        method = fetch_runner_logs
    elif handlername == 'fetch-codfreqs':
        method = fetch_codfreqs
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
