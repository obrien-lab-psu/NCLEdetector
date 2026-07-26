# Public Datastore R2 Release Plan

## Scope
Publish NCLEdetector_Datastore as a hybrid release:
- Exploded directory tree for per-file downloads
- Per-analysis archives next to each analysis directory
- Per-workflow archives next to each workflow directory
- Optional full-datastore archive

Release target:
- R2 bucket: ncle-public-data
- Public domain: https://downloads.ncleweb.org
- Catalog page: https://ncleweb.org/downloads
- Versioned prefix: NCLEdetector_Datastore/v1.0.0/

## Operational boundaries
- Repository work (plans, scripts, docs, validation helpers) is tracked and may require git commits in this repository.
- Publish payload work (release content to host) is based on an external filesystem tree, not the git working tree:
	- /scratch/ims86/NCLEdetector_Datastore
- We should treat these as two coupled tracks:
	- Track A: repository changes with normal review/commit checkpoints
	- Track B: datastore preparation and upload validation outside git

## Working mode
We execute one step at a time, verify output, then proceed.
No public exposure until private upload validation passes.

## Phase 1: Freeze and validate release payload

### Step 1. Freeze release identifiers
Decide and record:
- dataset version (example: v1.0.0)
- release date
- software version
- git commit hash
- datastore root path

Deliverable:
- release metadata values confirmed in writing

### Step 2. Finalize release tree shape
Under datastore root, confirm required top-level files and workflow layout.
Recommended root files:
- README.md
- LICENSE
- CITATION.cff
- RELEASE.json
- MANIFEST.tsv
- MANIFEST.json
- SHA256SUMS

Deliverable:
- agreed directory structure and archive placement

### Step 3. Pre-publication payload audit
Run and save outputs for:
- total size
- file count
- largest files
- symlink listing
- empty files
- temporary/artifact patterns

Gate:
- all suspicious files reviewed

### Step 4. Archive integrity validation
Validate all tar.gz files:
- gzip integrity test
- tar listing test

Gate:
- zero failed archives

### Step 5. Manifest generation
Generate MANIFEST.tsv and MANIFEST.json from finalized payload.
JSON should include at least path, size_bytes, and type.

Gate:
- manifest row/object count matches payload expectations

### Step 6. Checksum generation
Generate SHA256SUMS after all payload files are finalized.
Exclude MANIFEST.tsv, MANIFEST.json, SHA256SUMS from checksum inputs.

Gate:
- checksum file generated and spot-verified

## Phase 2: Private R2 staging

### Step 7. Create private bucket and token
Cloudflare:
- enable R2
- create bucket ncle-public-data (private)
- keep r2.dev disabled
- create bucket-scoped Object Read & Write token

Capture securely:
- Access Key ID
- Secret Access Key
- Account ID
- endpoint URL

Gate:
- credentials available and stored securely

### Step 8. Configure rclone remote
Configure remote (example name: r2-ncle) with:
- type s3
- provider Cloudflare
- region auto
- endpoint https://ACCOUNT_ID.r2.cloudflarestorage.com
- no_check_bucket true

Gate:
- rclone lsf remote bucket succeeds

### Step 9. Define versioned destination
Set:
- MASTER=/scratch/ims86/NCLEdetector_Datastore
- BUCKET=r2-ncle:ncle-public-data
- VERSION=v1.0.0
- DEST=${BUCKET}/NCLEdetector_Datastore/${VERSION}

Gate:
- source/destination printed and reviewed

### Step 10. Pilot upload and round-trip verify
Upload README.md to DEST and download back.
Compare local and downloaded files byte-for-byte.

Gate:
- round-trip compare succeeds

### Step 11. Full dry-run
Run rclone copy --dry-run with explicit excludes.
Inspect for path nesting and unexpected payload.

Gate:
- dry-run report approved

### Step 12. Full private upload
Run full rclone copy with conservative transfer settings.
Preserve upload log.

Gate:
- no critical errors in upload log

### Step 13. Private upload validation
Run:
- rclone size local vs remote
- rclone check --size-only
- remote tree inspection
- archive spot-download + checksum compare

Gate:
- all checks pass before public access

## Phase 3: Public access + catalog

### Step 14. Configure CORS
Allow origins:
- https://ncleweb.org
- https://www.ncleweb.org
Methods: GET, HEAD
Headers: Range
Expose: Content-Length, Content-Range, ETag, Last-Modified

Gate:
- CORS response headers verified with Origin request

### Step 15. Attach custom domain
Attach downloads.ncleweb.org to bucket and wait for Active.
Keep r2.dev disabled.

Gate:
- domain status Active

### Step 16. Public download tests
Validate:
- small file HEAD/GET
- archive HEAD
- range request
- CORS check on MANIFEST.json

Gate:
- all endpoint checks pass

### Step 17. Publish manifest-driven catalog
At https://ncleweb.org/downloads:
- load MANIFEST.json
- searchable file/path listing
- folder tree navigation
- mark .tar.gz as complete bundle downloads
- copy direct URL support
- show release metadata, citation, checksum instructions

Gate:
- browse/search/download UX verified

## Phase 4: Launch and update policy

### Step 18. Launch checklist
Confirm all technical gates passed and publish release announcement.

### Step 19. Future updates
Do not mutate v1.0.0 payload silently.
Publish new versions under new prefixes, example v1.1.0.
Rotate or revoke write token when release window closes.

## Immediate next actions
1. Collect Step 1 release identifiers from the current repository state.
2. Confirm the datastore publish root is /scratch/ims86/NCLEdetector_Datastore.
3. Decide whether full-datastore archive is included in v1.0.0.
4. Set a repository commit checkpoint before running release automation steps.

## Step 1 snapshot (captured 2026-07-26)
- release_date: 2026-07-26
- git_commit: 32096314b8940f7583b5f02ae817025bc2b4f1de
- git_branch: main
- software_version (from NCLEdetector/__init__.py): 1.1.6
- datastore_root: /scratch/ims86/NCLEdetector_Datastore
- datastore_size: 1.6T
- repository_state: dirty (34 changes currently present)

Checkpoint note:
- Because there are many in-progress repository edits, create a checkpoint commit (or a clearly labeled WIP commit sequence) before running any release automation that modifies scripts/configs.
