NIM = nim c

conduit:
	cd poaV2 && make
	$(NIM) --threads:on --passL:poaV2/liblpo.a --passL:poaV2/align_score.o conduit.nim 

