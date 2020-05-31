build:
	@docker build . -t hivdb/codfreq:latest-5prime

release: build
	@docker push hivdb/codfreq:latest-5prime

restart-ray: stop-ray start-ray

start-ray:
	@pipenv run ray start --head --redis-password 2g^jEK6O!

stop-ray:
	@pipenv run ray stop

debug:
	@docker run -it --rm hivdb/codfreq:latest-5prime bash
