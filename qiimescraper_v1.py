#!/usr/bin/env python3
    parser.add_argument("--r1-suffix", default="_R1_001.fastq.gz", help="Forward read filename suffix.")
    parser.add_argument("--r2-suffix", default="_R2_001.fastq.gz", help="Reverse read filename suffix.")
    parser.add_argument(
        "--sample-contains",
        help="Optional substring filter for sample FASTQ names. Useful for processing one marker at a time.",
    )

    parser.add_argument("--trunc-len-f", type=int, default=0, help="DADA2 forward truncation length.")
    parser.add_argument("--trunc-len-r", type=int, default=0, help="DADA2 reverse truncation length.")
    parser.add_argument("--trim-left-f", type=int, default=0, help="DADA2 trim-left forward bases.")
    parser.add_argument("--trim-left-r", type=int, default=0, help="DADA2 trim-left reverse bases.")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = Path(args.manifest) if args.manifest else outdir / "manifest.tsv"

    if args.full:
        run_full_pipeline(args)
        return

    if not args.step:
        parser.error("Choose either --full or --step")

    if args.step == "manifest":
        generate_manifest(
            fastq_dir=Path(args.fastq_dir),
            manifest_path=manifest,
            r1_suffix=args.r1_suffix,
            r2_suffix=args.r2_suffix,
            sample_contains=args.sample_contains,
        )
        return

    if args.step == "import":
        step_import(manifest, outdir)
        return

    if args.step == "cutadapt":
        forward_primers, reverse_primers = get_primers(args)
        step_cutadapt(outdir, forward_primers, reverse_primers, args.cores, args.discard_untrimmed)
        return

    if args.step == "dada2":
        step_dada2(
            outdir=outdir,
            cores=args.cores,
            trunc_len_f=args.trunc_len_f,
            trunc_len_r=args.trunc_len_r,
            trim_left_f=args.trim_left_f,
            trim_left_r=args.trim_left_r,
        )
        return

    if args.step == "taxonomy":
        step_taxonomy(outdir, Path(args.classifier), args.cores)
        return

    if args.step == "barplot":
        step_barplot(outdir, Path(args.metadata))
        return


if __name__ == "__main__":
    main()