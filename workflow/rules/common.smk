# Shared definitions: sample sheet, wildcard constraints, helper functions.
OUT = config["outdir"].rstrip("/")
SCRIPTS = os.path.join(workflow.basedir, "scripts")
ENVS = os.path.join(workflow.basedir, "envs")

SAMPLES = mb.load_samples(config["samples"], config["design"]["diagnosis_levels"])
SIDX = SAMPLES.set_index("sample_id")

_real = SIDX[SIDX.sample_type == "sample"]
MGX = _real[_real.data_type == "MGX"].index.tolist()
S16_REAL = _real[_real.data_type == "16S"].index.tolist()
S16_ALL = SIDX[SIDX.data_type == "16S"].index.tolist()          # includes negative controls
ASM = [s for s in MGX if SIDX.loc[s, "assemble"]]
BATCHES_16S = sorted(set(SIDX.loc[S16_ALL, "batch"]))

if not MGX or not S16_REAL:
    raise WorkflowError("Both 16S and MGX samples are required (the project compares the two platforms).")


def _alt(ids):
    return "(" + "|".join(re.escape(i) for i in ids) + ")" if ids else "(?!)"


wildcard_constraints:
    sample=r"[A-Za-z0-9_.\-]+",
    mate="[12]",
    layer=_alt(list(config["layers"])),
    method="(clr_lmm|maaslin2|ancombc2|aldex2)",
    clade=r"t__SGB[0-9A-Za-z_]+",


MGX_RE, S16_RE, ASM_RE = _alt(MGX), _alt(S16_ALL), _alt(ASM)
BATCH_RE = _alt(BATCHES_16S)
USE_HUMANN = config["shotgun"]["humann"]["enabled"]
USE_KRAKEN = config["shotgun"]["kraken2"]["enabled"]
USE_MAGS = config["mags"]["enabled"] and len(ASM) > 0
CLADES = list(config["shotgun"]["strainphlan"]["clades"])
PRIMERS = bool(config["amplicon"]["primer_fwd"] and config["amplicon"]["primer_rev"])

# layers that can be built with the configured options
ACTIVE_LAYERS = [l for l in config["layers"] if not (l == "mgx_pathway" and not USE_HUMANN)]
DA_JOBS = [(l, m) for l in ACTIVE_LAYERS for m in config["layers"][l]["da_methods"]]


def raw_reads(wc):
    return {"r1": f"{OUT}/raw/{wc.sample}_R1.fastq.gz", "r2": f"{OUT}/raw/{wc.sample}_R2.fastq.gz"}


def clean_reads(sample):
    return [f"{OUT}/shotgun/clean/{sample}_R1.fastq.gz", f"{OUT}/shotgun/clean/{sample}_R2.fastq.gz"]
