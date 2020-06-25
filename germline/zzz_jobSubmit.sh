#!/bin/bash

bamfiles=$1

bsub < bamSplice.bsub $bamfiles

