NIM = nim c

conduit:
	cd src/poaV2 && make
	mkdir bin/
	cd src/ && $(NIM) -d:release --threads:on --passL:poaV2/liblpo.a --passL:poaV2/align_score.o conduit.nim && mv src/conduit bin/
	cd src/ && $(NIM) -d:release conduitUtils.nim && mv src/conduitUtils bin/

