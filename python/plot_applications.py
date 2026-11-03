"""Plot the empirical and model spectra with pointwise reference envelopes."""

from pathlib import Path
import argparse
from manuscript_outputs import application_figure, reference_envelopes, style


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("figures/Fig4.pdf"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    reference_envelopes()
    style()
    application_figure(args.output.parent,args.output.name)


if __name__ == "__main__":
    main()
