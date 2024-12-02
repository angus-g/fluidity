#!/bin/bash

awk '/INTEGER\*8 :: /{print "       " $$0}' $1 >> $2
awk '/PARAMETER /{print "       " $$0}' $1 >> $2
awk '/INTEGER\*8 FUNCTION/{print "       " $$1 " " $$3}' $1 >> $2
