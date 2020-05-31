build:
	@docker build . -t hivdb/codfreq:latest-5prime

release: build
	@docker push hivdb/codfreq:latest-5prime

run-redis:
	@docker rm -f chiro-devredis 2>/dev/null || true
	@docker run \
		-d --name=chiro-devredis \
		--publish 127.0.0.1:16379:6379 \
		--volume=$(shell pwd)/redis-data:/data \
		redis:5 redis-server --appendonly yes

run-ray: run-redis
	@sleep 10
	@pipenv run ray --address 16379

debug:
	@docker run -it --rm hivdb/codfreq:latest-5prime bash
