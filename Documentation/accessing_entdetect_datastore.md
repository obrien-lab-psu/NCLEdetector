Below is a README-ready Markdown block you can give users. I’d recommend **GoCommands (`gocmd`)** as the primary route because CyVerse documents it as a lightweight cross-platform command-line tool for Data Store transfer/sync, including checksum verification, diff-based syncing, and multithreaded transfers. ([CyVerse Learning Materials][1])

````markdown
# Downloading the NCLEdetector CyVerse Data Store Dataset from the Command Line

This dataset is hosted on the CyVerse Data Store:

```text
/iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore
````

For large downloads, **do not use the browser download button**. Use the command line so the transfer can be resumed, verified, and run in a terminal multiplexer such as `tmux` or `screen`.

---

## Recommended method: CyVerse GoCommands

CyVerse GoCommands (`gocmd`) is the recommended command-line tool for downloading large folders from the CyVerse Data Store.

It supports:

* recursive folder downloads
* progress display
* checksum verification
* resuming/skipping files that already downloaded
* parallel transfer threads

---

## 1. Install GoCommands

### Option A: Linux or macOS install script

```bash
curl -fsSL https://raw.githubusercontent.com/cyverse/gocommands/main/install_gocmd.sh | bash
```

Restart your shell or make sure the installed `gocmd` binary is on your `PATH`.

Check that it works:

```bash
gocmd --help
```

### Option B: Conda install

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install gocommands
```

---

## 2. Configure access

### Option A: Use your CyVerse account

Run:

```bash
gocmd init
```

Use the following values when prompted:

```text
host: data.cyverse.org
port: 1247
username: <your CyVerse username>
zone: iplant
password: <your CyVerse password>
```

Test access:

```bash
gocmd ls /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore
```

### Option B: Anonymous/public access

If the dataset has been shared with anonymous read access, you can configure access using:

```text
username: anonymous
password: leave blank
zone: iplant
host: data.cyverse.org
port: 1247
```

Run:

```bash
gocmd init
```

Then test:

```bash
gocmd ls /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore
```

If anonymous access fails, use a CyVerse account or confirm that the folder has been shared with `anonymous` read permissions, not only the CyVerse `public` user.

---

## 3. Start the download

Choose a local destination with plenty of free space.

For example:

```bash
mkdir -p /path/to/local/downloads
```

Then run:

```bash
gocmd get \
  --progress \
  -k \
  --diff \
  --thread_num 8 \
  /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore \
  /path/to/local/downloads/
```

This will create:

```text
/path/to/local/downloads/NCLEdetector_Datastore/
```

### What the options mean

```text
--progress       Show transfer progress
-k               Verify file integrity using checksums
--diff           Skip files that already exist and only transfer new/changed files
--thread_num 8   Use 8 parallel transfer threads
```

For slow or unstable connections, try fewer threads:

```bash
--thread_num 4
```

For high-bandwidth institutional/HPC connections, you can try:

```bash
--thread_num 10
```

Avoid setting the thread count extremely high, because too many simultaneous transfers can overload your system or hit connection limits.

---

## 4. Resume an interrupted download

If the transfer fails or your connection drops, rerun the exact same command:

```bash
gocmd get \
  --progress \
  -k \
  --diff \
  --thread_num 8 \
  /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore \
  /path/to/local/downloads/
```

Because `--diff` is enabled, files that already downloaded successfully should be skipped, and only missing or changed files should be transferred.

---

## 5. Recommended: run inside `tmux` or `screen`

For very large downloads, run the command inside `tmux` so it continues even if your SSH session disconnects.

```bash
tmux new -s ncledetector_download
```

Run the `gocmd get ...` command inside the session.

Detach from the session with:

```text
Ctrl-b d
```

Reconnect later with:

```bash
tmux attach -t ncledetector_download
```

---

## 6. Optional: download to scratch space on an HPC system

On an HPC or cloud instance, download to a high-capacity scratch directory rather than your home directory.

Example:

```bash
mkdir -p /scratch/$USER/NCLEdetector_download

gocmd get \
  --progress \
  -k \
  --diff \
  --thread_num 8 \
  /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore \
  /scratch/$USER/NCLEdetector_download/
```

---

## 7. WebDAV fallback

If GoCommands is not available, CyVerse also exposes Data Store paths through WebDAV.

Base URL:

```text
https://data.cyverse.org/dav
```

Dataset path:

```text
https://data.cyverse.org/dav/iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore
```

For anonymous access, CyVerse WebDAV may be accessed using:

```text
username: anonymous
password: leave blank
```

However, for very large recursive downloads, `gocmd` is strongly preferred over WebDAV because it supports Data Store-aware sync, checksum verification, and parallel transfer options.

---

## Troubleshooting

### Permission denied

Try logging in with a CyVerse account instead of anonymous access.

Also confirm that the folder has been shared with the appropriate read permissions.

### Transfer interrupted

Rerun the same `gocmd get --diff ...` command. Existing files should be skipped.

### Download is slow

Try changing the thread count:

```bash
--thread_num 4
```

or

```bash
--thread_num 10
```

The best value depends on your network, storage speed, and whether you are downloading many small files or fewer large files.

### Checksums make the transfer slower

The `-k` flag verifies file integrity and is recommended for scientific datasets. If you only need a quick exploratory copy, you can omit `-k`, but for final analyses, checksum verification is recommended.

---

## Minimal command

For most users, this is the command to run:

```bash
gocmd get \
  --progress \
  -k \
  --diff \
  --thread_num 8 \
  /iplant/home/shared/NCEMS/NCLEdetector/NCLEdetector_Datastore \
  /path/to/local/downloads/
```

```

The main thing I’d emphasize to users is: **rerun the same `gocmd get --diff ...` command after failures**. CyVerse’s GoCommands docs describe `--diff` as transferring only new or modified files by comparing sizes/checksums, while `-k` enables checksum verification and `--thread_num` enables parallel transfers. :contentReference[oaicite:1]{index=1}
::contentReference[oaicite:2]{index=2}
```

[1]: https://learning.cyverse.org/ds/gocommands/ "Manage Your Data with GoCommands - CyVerse Learning Materials"

