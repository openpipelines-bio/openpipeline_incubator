#!/bin/bash

## VIASH START
## VIASH END

cat > _viash.yaml << EOM
info:
  test_resources:
    - type: s3
      path: s3://openpipelines-bio/openpipeline_incubator/resources_test
      dest: resources_test
EOM

echo ">> Run aws s3 sync"
"$meta_executable" \
  --input _viash.yaml \
  --output . \
  --quiet

echo ">> Check whether the right files were copied"
[ ! -f resources_test/qc_sample_data/sample_one.qc.h5mu ] && echo file should have been copied && exit 1

echo ">> Test succeeded!"
