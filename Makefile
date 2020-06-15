build: requirements.txt
	@docker build . -t hivdb/codfreq-runner:latest

requirements.txt: Pipfile Pipfile.lock
	@pipenv lock --requirements > requirements.txt

release: build
	@docker push hivdb/codfreq-runner:latest

syncrefs:
	@pipenv run aws s3 sync refs s3://codfreq-assets.hivdb.org/refs

debug:
	@docker run -it --rm \
		--volume ~/.aws:/root/.aws:ro \
		hivdb/codfreq-runner:latest bash
