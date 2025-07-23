SUBDIRS = obstools

.PHONY: all $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
