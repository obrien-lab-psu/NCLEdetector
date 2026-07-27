# Maintaining the NCLEdetector Public Datastore Release

This document describes how to update existing archives in the Cloudflare R2 bucket and add new analysis results to the public release.

---

## Overview

The NCLEdetector datastore is published as versioned archives to Cloudflare R2:

- **Bucket:** `ncle-public-data`
- **Public domain:** `https://downloads.ncleweb.org`
- **Versioned prefix:** `NCLEdetector_Datastore/v{VERSION}/`
- **Archives uploaded:** 18 per-analysis .tar.gz + 1 user_input.tar.gz + metadata files

Each version is immutable — do not modify files in place. Publish new versions under new prefixes (e.g., `v1.0.1`, `v1.1.0`).

---

## Adding New Archives to an Existing Release

### Scenario: Workflow2/OP.tar.gz finally completed, and you want to add it to v1.0.0

#### Step 1: Verify the new archive locally

```bash
ARCHIVE="/scratch/ims86/NCLEdetector_Datastore/outputs/workflow2/OP.tar.gz"

# Check size and completion
ls -lh "$ARCHIVE"

# Verify tar integrity
tar -tzf "$ARCHIVE" | head -20
echo "---"
tar -tzf "$ARCHIVE" | tail -20

# Compute checksum for records
sha256sum "$ARCHIVE"
```

#### Step 2: Upload to R2

```bash
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

ARCHIVE="/scratch/ims86/NCLEdetector_Datastore/outputs/workflow2/OP.tar.gz"
DEST="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0"

# Dry-run first
rclone copyto "$ARCHIVE" "$DEST/outputs/workflow2/OP.tar.gz" --dry-run

# Upload
rclone copyto "$ARCHIVE" "$DEST/outputs/workflow2/OP.tar.gz" --progress

# Verify round-trip
echo "Local checksum:"
md5sum "$ARCHIVE"

echo "Remote checksum:"
rclone cat "$DEST/outputs/workflow2/OP.tar.gz" | md5sum
```

#### Step 3: Update release metadata

If you added new archives to v1.0.0, regenerate manifests and checksums:

```bash
cd /scratch/ims86/NCLEdetector_Datastore

# Regenerate manifests (now includes OP.tar.gz)
python3 scripts/generate_release_manifests.py

# Regenerate checksums
sbatch assets/slurm/scripts/generate_release_checksums.slurm
# (or run locally: python3 scripts/generate_release_checksums.py)
```

#### Step 4: Upload updated metadata

```bash
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

DEST="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0"
MASTER="/scratch/ims86/NCLEdetector_Datastore"

# Upload updated manifests and checksums
rclone copyto "$MASTER/MANIFEST.tsv" "$DEST/MANIFEST.tsv" --progress
rclone copyto "$MASTER/MANIFEST.json" "$DEST/MANIFEST.json" --progress
rclone copyto "$MASTER/SHA256SUMS" "$DEST/SHA256SUMS" --progress

# Verify sizes
rclone size "$DEST"
du -sh "$MASTER" --exclude='outputs/workflow2/OP.tar.gz'  # Old payload size
du -sh "$MASTER"  # New total
```

---

## Publishing a New Release Version

### Scenario: You have completed v1.1.0 with new analysis results

#### Step 1: Prepare release metadata

```bash
DATASTORE="/scratch/ims86/NCLEdetector_Datastore"
VERSION="v1.1.0"

# Update or create new release files
cat > "$DATASTORE/RELEASE.json" << 'EOF'
{
  "version": "1.1.0",
  "release_date": "2026-07-27",
  "software_version": "1.2.0",
  "git_commit": "abc123def456",
  "description": "Updated release with new analyses",
  "download_url": "https://downloads.ncleweb.org/NCLEdetector_Datastore/v1.1.0/",
  "archive_count": 18,
  "total_size_bytes": 1589000000000,
  "workflows": ["workflow1", "workflow2", "workflow3", "workflow4"],
  "license": "GNU GPL v3"
}
EOF

# Update README if needed
# Update CITATION.cff with new version
```

#### Step 2: Generate manifests and checksums for v1.1.0

```bash
cd "$DATASTORE"

# Create separate directory for v1.1.0 prep if desired, or regenerate in place
python3 scripts/generate_release_manifests.py
sbatch assets/slurm/scripts/generate_release_checksums.slurm
```

#### Step 3: Create and submit upload SLURM job for v1.1.0

Create a new upload script or modify existing:

```bash
# Copy the existing upload script and modify destination
cp assets/slurm/scripts/upload_release_to_r2.slurm \
   assets/slurm/scripts/upload_release_v1.1.0_to_r2.slurm

# Edit the new script to change VERSION variable:
# DEST="${BUCKET}:${BUCKET_NAME}/NCLEdetector_Datastore/v1.1.0"
# (and update any archive exclusion if needed)

# Submit
sbatch assets/slurm/scripts/upload_release_v1.1.0_to_r2.slurm
```

#### Step 4: Validate remote v1.1.0 payload

```bash
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

DEST_NEW="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.1.0"

# Size check
rclone size "$DEST_NEW"

# Manifest check
rclone cat "$DEST_NEW/MANIFEST.json" | python3 -m json.tool | head -20

# Sample file integrity
curl -I "https://downloads.ncleweb.org/NCLEdetector_Datastore/v1.1.0/README.md"
```

#### Step 5: Update download instructions

Edit [index.md](index.md) `### 2. Fetch and set the datastore path` section to point users to the new version URL.

#### Step 6: Archive v1.0.0 in release notes

Document the v1.0.0 release in CHANGELOG.md or release notes for future reference.

---

## Updating Individual Archives (Advanced)

### Keep v1.0.0 immutable; replace only on major updates

**Do NOT:** directly delete and re-upload a single archive under the same name/version.

**Instead:** Create a new version (v1.0.1) and re-upload the corrected archive.

```bash
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

# Old version (immutable — do not modify)
# DEST_OLD="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0"

# New version (for bug fixes or updates)
DEST_NEW="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.1"

# Upload corrected archive
rclone copyto /path/to/corrected_archive.tar.gz "$DEST_NEW/outputs/workflow2/corrected_archive.tar.gz"

# Regenerate and re-upload manifests/checksums for v1.0.1
python3 scripts/generate_release_manifests.py
sbatch assets/slurm/scripts/generate_release_checksums.slurm
rclone copyto MANIFEST.tsv "$DEST_NEW/MANIFEST.tsv"
rclone copyto MANIFEST.json "$DEST_NEW/MANIFEST.json"
rclone copyto SHA256SUMS "$DEST_NEW/SHA256SUMS"
```

---

## Revoking Upload Tokens After Release

Once a release version is finalized and no further uploads are expected for that version:

1. **In Cloudflare Dashboard:** Go to **R2 → API Tokens**
2. **Revoke or rotate** the bucket-scoped token used for uploads
3. **Document** in release notes or changelog that this token is retired

This prevents accidental modifications to published versions.

---

## Verification Checklist

After any upload:

- [ ] Remote size matches expected total
- [ ] Sample files download successfully (200 OK)
- [ ] Archive checksums match local→remote round-trip
- [ ] MANIFEST.json is valid JSON and loadable
- [ ] CORS headers present for cross-origin requests
- [ ] Public domain (downloads.ncleweb.org) serves files correctly
- [ ] Version prefix in URL is correct

---

## Support & Questions

- **Manifest format:** See [MANIFEST.json structure](MANIFEST.md) for schema details
- **Archive integrity:** Use `sha256sum -c SHA256SUMS` to verify downloads
- **Upload logs:** Check `/tmp/r2_upload_*.log` for detailed transfer information
- **rclone config:** Stored in `~/.config/rclone/rclone.conf` (keep secure!)

