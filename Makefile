requirements.txt: Pipfile Pipfile.lock
	@pipenv lock --requirements > requirements.txt
	# @sed -i '' -E 's/; (sys_platform|python_(full_)?version).*$$//' requirements.txt

login-ecr:
	@aws ecr get-login-password --region us-west-2 | docker login \
		--username AWS \
		--password-stdin 931437602538.dkr.ecr.us-west-2.amazonaws.com

build-runner: requirements.txt
	@docker build . -t hivdb/codfreq-runner:latest

release-runner: build-runner login-ecr deploy-profiles
	@docker push hivdb/codfreq-runner:latest
	@docker tag hivdb/codfreq-runner:latest 931437602538.dkr.ecr.us-west-2.amazonaws.com/hivdb/codfreq-runner:latest
	@docker push 931437602538.dkr.ecr.us-west-2.amazonaws.com/hivdb/codfreq-runner:latest

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

debug-runner:
	@docker run -it --rm \
		--volume $(PWD)/local:/local:rw \
		--volume ~/.aws:/root/.aws:ro \
		hivdb/codfreq-runner:latest bash

deploy-profiles:
	@aws s3 sync profiles s3://codfreq-assets.hivdb.org/profiles --delete

.PHONY: login-ecr *-runner *-controller deploy-profiles
