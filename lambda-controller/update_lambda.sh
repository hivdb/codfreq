#! /bin/bash

set -e

cd $(dirname $0)

. config.sh

aws lambda update-function-code \
  --region ${AWS_REGION} \
  --function-name ${FUNCTION_PREFIX}Dispatch \
  --publish \
  --zip-file ${1}
