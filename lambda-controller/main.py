# type: ignore
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


def save_file(uniqkey, config_file, content_type, payload):
    S3.put_object(
        Bucket=S3_BUCKET,
        Key='tasks/{}/{}'.format(uniqkey, config_file),
        ContentType=content_type,
        Body=payload.encode('U8')
    )


def save_pairinfo(uniqkey, pairinfo):
    save_file(
        uniqkey,
        'pairinfo.json',
        'application/json',
        json.dumps(pairinfo)
    )


def verify_pairinfo(pairinfo, fastq_files):
    known_names = set([])
    for idx, one in enumerate(pairinfo):
        if 'name' not in one:
            return False, "item {} misses property 'name'".format(idx)
        if 'pair' not in one:
            return False, "item {} misses property 'pair'".format(idx)
        if 'n' not in one:
            return False, "item {} misses property 'n'".format(idx)
        if one['name'] in known_names:
            return False, "item {} has duplicate name".format(idx)
        if len(one['pair']) > 2:
            return (
                False,
                "item {}'s property 'pair' has too many files".format(idx)
            )
        if len(one['pair']) < 2:
            return (
                False,
                "item {}'s property 'pair' has too few file".format(idx)
            )
        if one['n'] == 2 and any(f is None for f in one['pair']):
            return (
                False,
                "item {}'s property 'pair' doesn't match n=2".format(idx)
            )
        if (
            one['n'] == 1 and
            (one['pair'][0] is None or
             one['pair'][1] is not None)
        ):
            return (
                False,
                "item {}'s property 'pair' doesn't match n=1".format(idx)
            )
        files = {f for f in one['pair'] if f is not None}
        missing_files = files - fastq_files
        if missing_files:
            return (
                False,
                "item {}'s property 'pair' lists following files "
                "which are not uploaded: {}"
                .format(idx, ', '.join(missing_files))
            )
    return True


def load_taskmeta(uniqkey):
    obj = S3.get_object(
        Bucket=S3_BUCKET,
        Key='tasks/{}/taskmeta.json'.format(uniqkey)
    )
    data = json.load(obj['Body'])
    data['lastUpdatedAt'] = parse_dt(data['lastUpdatedAt'])
    return data


def verify_profiles(profiles):
    for profile in profiles:
        try:
            S3.head_object(
                Bucket=S3_BUCKET,
                Key='profiles/{}'.format(profile)
            )
        except S3.exceptions.NoSuchKey:
            return False
    return True


def get_or_default(map, key, default, type_):
    value = map.get(key)
    if value is None:
        return default
    return type_(value)


def save_fastp_config(uniqkey, opts):
    payload = json.dumps({
        'include_unmerged': get_or_default(
            opts, 'includeUnmerged', True, bool),
        'qualified_quality_phred': get_or_default(
            opts, 'qualifiedQualityPhred', 15, int),
        'unqualified_percent_limit': get_or_default(
            opts, 'unqualifiedPercentLimit', 40, int),
        'n_base_limit': get_or_default(
            opts, 'nBaseLimit', 5, int),
        'average_qual': get_or_default(
            opts, 'averageQual', 0, int),
        'length_required': get_or_default(
            opts, 'lengthRequired', 15, int),
        'length_limit': get_or_default(
            opts, 'lengthLimit', 0, int),
        'adapter_sequence': get_or_default(
            opts, 'adapterSequence', 'auto', str),
        'adapter_sequence_r2': get_or_default(
            opts, 'adapterSequenceR2', 'auto', str),
        'disable_adapter_trimming': get_or_default(
            opts, 'disableAdapterTrimming', False, bool),
        'disable_trim_poly_g': get_or_default(
            opts, 'disableTrimPolyG', False, bool),
        'disable_quality_filtering': get_or_default(
            opts, 'disableQualityFiltering', False, bool),
        'disable_length_filtering': get_or_default(
            opts, 'disableLengthFiltering', False, bool)
    })
    save_file(uniqkey, 'fastp-config.json', 'application/json', payload)


def save_cutadapt_config(uniqkey, opts):
    config = json.dumps({
        'error_rate': get_or_default(
            opts, 'errorRate', 0.1, float),
        'no_indels': get_or_default(
            opts, 'noIndels', True, bool),
        'times': get_or_default(
            opts, 'times', 1, int),
        'min_overlap': get_or_default(
            opts, 'minOverlap', 3, int)
    })
    three_end = '\n'.join(
        '>{}\n{}\n'.format(s['header'], s['sequence'])
        for s in opts['primerSeqs']
        if s['type'] == 'three-end'
    )
    five_end = '\n'.join(
        '>{}\n{}\n'.format(s['header'], s['sequence'])
        for s in opts['primerSeqs']
        if s['type'] == 'five-end'
    )
    both_end = '\n'.join(
        '>{}\n{}\n'.format(s['header'], s['sequence'])
        for s in opts['primerSeqs']
        if s['type'] == 'both-end'
    )
    save_file(uniqkey, 'cutadapt-config.json', 'application/json', config)
    if three_end:
        save_file(uniqkey, 'primers3.fa', 'text/x-fasta', three_end)
    if five_end:
        save_file(uniqkey, 'primers5.fa', 'text/x-fasta', five_end)
    if both_end:
        save_file(uniqkey, 'primers53.fa', 'text/x-fasta', both_end)


def save_ivar_config(uniqkey, opts):
    config = json.dumps({
        'min_length': get_or_default(opts, 'minLength', 0, int),
        'min_quality': get_or_default(opts, 'minQuality', 0, int),
        'sliding_window_width': get_or_default(
            opts, 'slidingWindowWidth', None, int),
        'include_reads_with_no_primers': get_or_default(
            opts, 'includeReadsWithNoPrimers', True, bool)
    })
    bedfp = StringIO()
    writer = csv.writer(bedfp, delimiter='\t')
    for row in opts.get('primerBeds') or []:
        start = get_or_default(row, 'start', None, int)
        end = get_or_default(row, 'end', None, int)
        name = get_or_default(row, 'name', None, str)
        strand = get_or_default(row, 'strand', None, str)
        if start is None or \
                end is None or \
                name is None or \
                strand not in ('+', '-'):
            # skip invalid rows
            continue
        writer.writerow([
            get_or_default(row, 'region', '<Unknown>', str),
            start,
            end,
            name,
            get_or_default(row, 'score', 60, int),
            strand
        ])
    save_file(uniqkey, 'ivar-trim-config.json', 'application/json', config)
    save_file(uniqkey, 'primers.bed', 'text/x-bed', bedfp.getvalue())


def create_task(request):
    """Step 1: create task meta file"""
    uniqkey = str(uuid.uuid4())
    now = save_taskmeta(uniqkey, 'created')
    options = request.get('options') if request else None
    if options:
        if 'fastpConfig' in options:
            save_fastp_config(uniqkey, options['fastpConfig'])
        if 'cutadaptConfig' in options:
            save_cutadapt_config(uniqkey, options['cutadaptConfig'])
        if 'ivarConfig' in options:
            save_ivar_config(uniqkey, options['ivarConfig'])

    return {
        'taskKey': uniqkey,
        'lastUpdatedAt': now,
        'status': 'created'
    }, 200


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
    location = S3.generate_presigned_url(
        'get_object',
        Params={
            'Bucket': S3_BUCKET,
            'Key': 'tasks/{}/response.json'.format(uniqkey)
        },
        ExpiresIn=300)

    return {
        'location': location
    }, 302


def fetch_codfreqs_zip(request):
    """Step 5 alt: fetch codfreq zip file"""
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
    location = S3.generate_presigned_url(
        'get_object',
        Params={
            'Bucket': S3_BUCKET,
            'Key': 'tasks/{}/codfreqs.zip'.format(uniqkey)
        },
        ExpiresIn=300)

    return {
        'location': location
    }, 302


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
            logStreamName='ecs/{}'.format(task_id),
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
    pairinfo = request.get('pairInfo')
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
        MaxKeys=len(fastq_files) + 6,
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

    all_profiles = set()
    for runner in runners:
        if not isinstance(runner, dict):
            return {'error': "invalid format of 'runners'"}, 400
        profile = runner.get('profile')
        all_profiles.add(profile)

    if not verify_profiles(all_profiles):
        return {
            'error': "invalid format of 'runners': unsupport 'profile'"
        }

    if pairinfo:
        if not verify_pairinfo(pairinfo, fastq_files):
            return {
                'error': "invalid format of 'pairInfo'"
            }
        save_pairinfo(uniqkey, pairinfo)

    task_ids = []
    for runner in runners:
        profile = runner['profile']
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
                        'bin/align-all-s3',
                        '-r', profile,
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
    elif handlername == 'fetch-codfreqs-zip':
        method = fetch_codfreqs_zip
    else:
        method = method_not_found
    response_body, status = method(request)

    common_headers = {
        'Access-Control-Allow-Origin': '*'
    }

    if status in (301, 302):
        return {
            'statusCode': status,
            'headers': {
                **common_headers,
                **response_body
            }
        }
    else:
        return {
            'statusCode': status,
            'headers': {
                **common_headers,
                'Content-Type': 'application/json'
            },
            'body': json.dumps(response_body)
        }
