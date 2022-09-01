# the user add the local file that will make up the library
SOURCES = reader reader_utils

CSRC = $(patsubst %,src/%.c,$(SOURCES))
COBJ = $(patsubst %,$(OBJDIR)/%.o,$(SOURCES))

all: preamble $(COBJ)

preamble:
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: src/%.c
	@echo compiling $<
	$(CC) $(CFLAGS) $(CINC) $(DEPINC) -c $< -o $@

clean:
	~rm -vf $(OBJDIR)
