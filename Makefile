build:
	@docker build . -t hivdb/codfreq:latest

release: build
	@docker push hivdb/codfreq:latest

debug:
	@docker run -it --rm hivdb/codfreq:latest bash
