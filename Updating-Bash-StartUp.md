# Updating Bash Startup — **“Explain‑It‑Like‑I’m‑Five” Edition** 

> **Goal:** Tell your machine where your bio‑toolbox lives so you can type
> `nanoDengue.sh` in *any* folder and it “just works.”

---

## 1. What *is* a “bash startup file”?

* **bash** = the program that runs in the background when you open a terminal
* **\~/.bashrc** = a little script bash reads **every** time it starts
* Whatever you put in `~/.bashrc` becomes the *default* behaviour for every new terminal window

Think of it as your coffee‑order sticky note for bash.

---

## 2. Two tiny lines that need to be added

```bash
export PATH="$HOME/bio-toolbox/bin:$PATH"
export BIO_CONF_DIR="$HOME/bio-configs"
```

**Plain‑English translation:**

| Line                    | Says to bash                                                                                             | Why it's needed                                                                                            |
| ----------------------- | -------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `export PATH=…`         | “Look for programs in the folder `~/bio-toolbox/bin` **first**, then everywhere else you normally look.” | So `nanoDengue.sh` (and future scripts) are found no matter where you are on the filesystem.              |
| `export BIO_CONF_DIR=…` | “Our global config files live in `~/bio-configs`.”                                                       | The script reads this variable to locate `nanoDengue.conf` (paths to minimap2, samtools, etc.). |

---

## 3. Step‑by‑step: copy–paste guide (with “where’s my .bashrc?”)

1. **Find & open your `.bashrc`:**

   * Location: it lives in your home folder as a *hidden* file called `~/.bashrc` (the `~` is shorthand for `/home/your‑username`).
   * First, check if it exists by entering `ls -a` *OR* `ls -a ~` (basically the same thing).
   * Command examples (one or more of these will already be on your machine):

     ```bash
     nano ~/.bashrc          # super‑simple editor (recommended for beginners. I opened it once and hated it though)
     vim  ~/.bashrc          # if you’re a vim person. might hurt your eyes if your .vimrc file isn't set-up for ease of vision
     emacs ~/.bashrc         # allows copy-paste on my machine, but might be be difficult to navigate
     code ~/.bashrc          # VS Code, if the `code` CLI is installed
     ```

     *(If the file doesn’t exist, the editor will start a blank file and create it when you save)*
2. **Scroll to the very bottom.**  Ignore everything above.
3. **Paste** the two‑line block from [Section 2](#2-two-tiny-lines-that-need-to-be-added) *exactly* as shown.
4. **Save & exit** the editor. See collapsible cheat sheet below!
    <details>
    <summary>Cheat-sheet: how to save & close in each editor</summary>

    | Editor | Save changes | Quit / Close |
    |--------|--------------|--------------|
    | **nano** | `Ctrl-O` then `Enter` | `Ctrl-X` |
    | **vim**  | `Shift-Z Z` *(saves **and** quits in one step w/o having to exit Insert Mode)* OR `:w` (after returning to Normal Mode) | ``:q` (or `:wq` to save *and* quit in one go if in Normal Mode) |
    | **emacs** | `Ctrl-X` `Ctrl-S` | `Ctrl-X` `Ctrl-C` |
    | **VS Code** (`code`) | `Ctrl-S` | `Ctrl-Q` *(Windows/Linux)* or close the tab/window |
    </details>

5. **Reload** the file so your current terminal picks up the changes:

   ```bash
   source ~/.bashrc
   ```

   *(Opening a new terminal tab does the same thing)*

---

## 4. Test That It Works!

```bash
which nanoDengue.sh   # should print /home/you/bio-toolbox/bin/nanoDengue.sh
printenv BIO_CONF_DIR # should print /home/you/bio-configs
```

If both commands show the expected paths, you’re done! Yay!

(If it's not working, *TELL ME!*)

---

## 5. FAQs (a.k.a. “But what if…?”)

| Question                                         | Quick answer                                                                                                                                                                  |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| “I get `command not found` for `nanoDengue.sh`.” | Double‑check that `nanoDengue.sh` really is inside `~/bio-toolbox/bin` *AND* that you copied the `export PATH=…` line correctly.                                              |
| “I have zsh, not bash.”                          | Edit `~/.zshrc` instead of `~/.bashrc`, then run `source ~/.zshrc`. Same two lines. Shouldn't be a problem since `.zshrc` is a MacOS specific file type.                                                                                          |
| “My configs live somewhere else.”                | Change the right‑hand path in `BIO_CONF_DIR` to match.                                                                                                                        |
| “I only want this for one project.”              | Skip the bashrc edits and, inside the project, run:<br>`export BIO_CONF_DIR=/path/to/configs`<br>`export PATH=/path/to/bio-toolbox/bin:$PATH`<br>before launching the script. |

That’s all there is to it—two lines, one reload, done!
