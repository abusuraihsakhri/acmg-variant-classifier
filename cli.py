#!/usr/bin/env python3
"""Command Line Interface for ACMG/AMP Variant Classifier & Functional Genomics.

Provides unified entry points for interactive mode, single-variant evaluation,
TSV/CSV/VCF batch processing, splice impact prediction, and pharmacogenomics (CPIC).
"""

import sys
from acmg_classifier.cli import (
    main,
    build_parser,
    _run_main_logic,
    _split_codes,
    _validate_input_path,
    _validate_output_path,
    _probability_float,
    run_interactive,
    run_batch_table,
    run_vcf,
)

if __name__ == "__main__":
    sys.exit(main())
