#!/usr/bin/env bash
set -e
export PATH=~/.local/bin:$PATH

source setup.sh
make setup
ls -l build

if [ "x${CI_RUN_TEST}" = 'xUNIT' ]
then
  make test
fi

if [ "x${CI_RUN_TEST}" = 'xINTEGRATION' ]
then
  ./.ci/integration-test.sh
fi
set +e
