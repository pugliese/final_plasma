#!/bin/bash

for k in {100..200}
do
  ./md.e p $k
  echo $k/200
done
