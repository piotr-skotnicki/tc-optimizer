#!/bin/sh

(cd isl && ./autogen.sh)

(cd pet && ./autogen.sh)

(cd barvinok && ./autogen.sh)

autoreconf -i
