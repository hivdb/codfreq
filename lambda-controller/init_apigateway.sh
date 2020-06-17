#! /bin/bash

set -e

cd $(dirname $0)

. config.sh

sed "s/\$\$FUNCPREFIX/${FUNCTION_PREFIX}/g;s/\$\$REGION/${AWS_REGION}/g" apigateway.json > /tmp/${FUNCTION_PREFIX}apigateway.json
trap "rm /tmp/${FUNCTION_PREFIX}apigateway.json" EXIT
aws apigateway import-rest-api \
  --region $AWS_REGION \
  --body "file:///tmp/${FUNCTION_PREFIX}apigateway.json"
