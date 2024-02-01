`nextflow run nf-core/pangenome -r 1.0.0 -profile docker -work-dir "path/to/work" -params-file "path/to/nf-params-Pan.yaml"`

`vg autoindex --prefix Pan-exome --workflow giraffe --gfa Pan-exome.fa.gz.gfaffix.unchop.Ygs.view.gfa`

`vg paths -S "REF" --list -x Pan-exome.giraffe.gbz > Pan-exome.reference_path.txt`
