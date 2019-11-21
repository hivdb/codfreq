build:
	@docker build . -t hivdb/codfreq:latest-5prime

release: build
	@docker push hivdb/codfreq:latest-5prime

debug:
	@docker run -it --rm hivdb/codfreq:latest-5prime bash
