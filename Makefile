build:
	@docker build . -t hivdb/codfish:latest

debug:
	@docker run -it --rm hivdb/codfish:latest bash
