#! /bin/bash

cd $(dirname $0)

. config.sh

LAMBDA_ARGS=$(get_lambda_args ${FUNCTION_PREFIX} ${ROLE} $1)

set -e

aws lambda create-function \
  --runtime python3.8 \
  --tags Category=${FUNCTION_PREFIX} \
  --role ${ROLE} \
  --publish \
  --zip-file ${1} \
  --region ${AWS_REGION} \
  --function-name ${FUNCTION_PREFIX}Dispatch \
  --handler main.dispatch \
  --description "Controller enterpoint of CodFreq Runner." \
  --timeout 20 \
  --memory-size 128
