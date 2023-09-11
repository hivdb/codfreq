get_lambda_args() {
    echo "--runtime python3.8 \
        --tags Category=${1} \
        --role ${2} \
        --publish \
        --zip-file ${3}"
}
