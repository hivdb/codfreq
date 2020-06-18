requirements.txt: Pipfile Pipfile.lock
	@pipenv lock --requirements > requirements.txt

build-runner: requirements.txt
	@docker build . -t hivdb/codfreq-runner:latest

release-runner: build-runner
	@docker push hivdb/codfreq-runner:latest

build-controller:
	@mkdir -p build
	@rm build/codfreq-controller.zip 2>/dev/null || true
	@zip -j build/codfreq-controller.zip lambda-controller/main.py

init-controller: build-controller
	@lambda-controller/init_lambda.sh fileb://$(realpath build/codfreq-controller.zip)

deploy-controller: build-controller
	@lambda-controller/update_lambda.sh fileb://$(realpath build/codfreq-controller.zip)

syncrefs:
	@pipenv run aws s3 sync refs s3://codfreq-assets.hivdb.org/refs

debug-docker:
	@docker run -it --rm \
		--volume ~/.aws:/root/.aws:ro \
		hivdb/codfreq-runner:latest bash


