# Makefile — build + versioned, structured install for comsa-rs.
#
# Layout produced by `make install`:
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/comsa-rs     the compiled binary
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/comsa        symlink to it (the pre-rename name)
#   $(PREFIX)/$(NAME)-$(VERSION)/README.md        documentation (see DOCS)
#   $(PREFIX)/$(NAME)-$(VERSION)/CHANGELOG.md
#   $(PREFIX)/$(NAME)-$(VERSION)/LICENSE
#
# Examples:
#   make                          # build (release)
#   make test                     # unit + integration tests
#   make install PREFIX=/opt      # -> /opt/comsa-rs-0.1.0/...
#   make dist PREFIX=/tmp/stage   # build a versioned tarball
#   make version                  # print the version from Cargo.toml

PREFIX     ?= /usr/local
DESTDIR    ?=
NAME       ?= comsa-rs
CARGO      ?= cargo
INSTALL    ?= install
BIN         = comsa-rs
# Installed as a symlink to $(BIN) so scripts that call the old name keep working.
ALIAS       = comsa
MANIFEST    = Cargo.toml
# Documentation shipped in a binary release. Listed by name rather than globbed as
# *.md, so that a maintainer note at the top level does not end up in a release.
DOCS        = README.md CHANGELOG.md LICENSE

# The version of record is [package] version in Cargo.toml; the sed range confines
# the match to that table.
VERSION := $(shell sed -n '/^\[package\]/,/^\[/s/^version *= *"\([^"]*\)".*/\1/p' $(MANIFEST) | head -1)
$(if $(VERSION),,$(error could not read version from $(MANIFEST)))

# Versioned package root (DESTDIR-aware for staged packaging).
PKGDIR = $(DESTDIR)$(PREFIX)/$(NAME)-$(VERSION)

.PHONY: all build test version install uninstall dist clean

all: build

version:
	@echo $(VERSION)

build:
	$(CARGO) build --release

test:
	$(CARGO) test --release

install: build
	$(INSTALL) -d $(PKGDIR)/bin
	$(INSTALL) -m 0755 target/release/$(BIN) $(PKGDIR)/bin/$(BIN)
	ln -sf $(BIN) $(PKGDIR)/bin/$(ALIAS)
	@for d in $(DOCS); do \
	  $(INSTALL) -m 0644 $$d $(PKGDIR)/ && echo "  doc       <- $$d"; \
	done
	@echo "Installed $(NAME) $(VERSION) -> $(PKGDIR)"

uninstall:
	rm -rf $(PKGDIR)

# Versioned tarball of the installed hierarchy: <name>-<version>.tar.gz
dist: install
	tar -C $(DESTDIR)$(PREFIX) -czf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	@echo "Wrote $(NAME)-$(VERSION).tar.gz"

clean:
	$(CARGO) clean
