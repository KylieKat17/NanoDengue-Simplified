# NanoDengue‑Simplified

*A one‑command wrapper around the original* **[NanoDengue](https://github.com/rajithadp/NanoDengue)** *pipeline—designed for a single wet-lab biologist who’d rather not fiddle with shell paths every run*

---

## ✨ What This Wrapper Adds

| Feature                         | Why you care                                                                                                                                      |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Global Config**               | Puts absolute paths to **minimap2**, **samtools**, **NanoPlot** in **one** file (`$HOME/bio-configs/nanoDengue.conf`). No more hard‑coding. |
| **Auto‑detect inputs**          | Drop a `fastq_pass/` folder and a `*.fasta` reference into any run directory, then just type `nanoDengue.sh` and go.                                  |
| **CLI Override Flags**               | Need a different reference or output folder? Use `-r`, `-i`, and `-o` flags to point at a different reference, reads folder, or NanoPlot output dir *without editing the script*              |
| **Inline Help & Sanity Checks** | The script tells you if a tool is missing *before* running and wasting your time, and prints a summary of which paths it will use.                                              |
| **Heavily‑Commented Source**    | Every block includes “why” comments, so making adjustments or debugging is (hopefully) more straightforward on a diff machine.                                                                    |

---

## 🚀 Quick Start (no Git CLI required)

1. **Download the repo/ZIP file**
   
   On the private GitHub page, click **`<> Code ➜ Download ZIP`**. Unzip anywhere (e.g. `~/Downloads`). If downloading via email, download and just unzip the same way!
   
3. **Move the toolbox + config into place**

   ```bash
   mv ~/Downloads/NanoDengueSimplified/bio-toolbox  $HOME/
   mv ~/Downloads/NanoDengueSimplified/bio-configs  $HOME/
   ```

   *Nothing to rename:* the repo already ships with `bio-configs/nanoDengue.conf` — open it and edit the right‑hand paths *only if* your tools aren’t on `$PATH`.
4. **Update your shell startup**
   
   Follow **[`Updating-Bash-StartUp`](./Updating-Bash-StartUp.md)** (two copy‑paste lines in `~/.bashrc`, then `source ~/.bashrc`).
6. **Prepare a run folder**

   ```text
   my‑run/
   ├─ fastq_pass/            # barcoded sub‑dirs / *.fastq.gz
   └─ dengue.fa             # reference genome
   ```

   *This is an example!* No clue how your file management is set up so this might differ!
7. **Run the pipeline**

   ```bash
   cd my-run
   nanoDengue.sh            # auto‑detects everything
   # Examples:
   nanoDengue.sh -r ZIKV.fa                # custom reference
   nanoDengue.sh -i reads -o QC_Plots      # custom reads & output folder
   ```

> **Tip:** If you *do* have Git installed, `git clone` works too. The folder names are identical, so Steps 2‑5 are unchanged.

---

## 🗄️  Toolbox & Config Layout

```text
$HOME/
├─ bio-toolbox/
│   ├─ bin/
│   │   └─ nanoDengue.sh      # symlink → ../pipelines/nanoDengue.sh
│   └─ pipelines/
│       └─ nanoDengue.sh      # the actual, commented script
└─ bio-configs/
    └─ nanoDengue.conf        # MINIMAP2=… SAMTOOLS=… NANOPLOT=…
```

### nanoDengue.conf template

```bash
# Only fill in a path if the tool is NOT already on $PATH.
MINIMAP2=/opt/minimap2/minimap2
SAMTOOLS=/usr/local/bin/samtools
NANOPLOT=$HOME/miniconda3/envs/nano/bin/NanoPlot
# GZIP=/usr/bin/gzip
# GUNZIP=/usr/bin/gunzip
```
> This is located in the `bio-configs` directory already!

---

## 🔬 What the script does (internals)

1. **NanoPlot QC**  → Kde/Hex/Dot plots for every `*.fastq.gz`.
2. **Decompress**   → writes `*.fastq` copies (keeps `.gz`).
3. **Alignment**    → `minimap2 -a reference reads > SAM`.
4. **BAM workflow** → sort → index → consensus FASTA (`samtools`).
5. **Outputs** per read file:

   * `Nano/<subdir>/<sample>/*.png|.html`  – QC plots
   * `<sample>_minimap2.sam`
   * `<sample>_minimap2_sorted.bam` + `.bai`
   * `<sample>_minimap2_sorted.consensus.fasta`

---

## 🧩 Troubleshooting Cheatsheet

| Symptom                        | Likely Fix                                                                                    |
| ------------------------------ | -------------------------------------------------------------------------------------- |
| `command not found minimap2`   | Add `MINIMAP2=/abs/path/minimap2` in `nanoDengue.conf` or install minimap2 on `$PATH`. |
| `Permission denied` on script  | `chmod +x $HOME/bio-toolbox/pipelines/nanoDengue.sh`                                   |
| Symlink issues on immutable FS | Run with `-i fastq_pass -o Nano` so symlinks aren’t needed.                            |
| Wrong reference auto‑picked    | Use `-r correctRef.fa`.                                                                |

More context is baked into the script’s comments.

---

## 📜 License & credits

* Original algorithm – **NanoDengue** by [rajithadp](https://github.com/rajithadp/NanoDengue) (MIT).
* Wrapper & docs – *created for a single private user*

