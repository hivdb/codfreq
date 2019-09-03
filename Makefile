build:
	@docker build . -t hivdb/codfish:latest

release: build
	@docker push hivdb/codfish:latest

debug:
	@docker run -it --rm hivdb/codfish:latest bash
