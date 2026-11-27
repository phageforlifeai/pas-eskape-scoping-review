#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# push_to_your_repo.sh
#
# Pushes the prepared repository to a GitHub repo you have already created.
#
# Usage:
#   bash push_to_your_repo.sh <repo-url> [token]
#
# Examples:
#   bash push_to_your_repo.sh https://github.com/yourname/pas-eskape-scoping-review.git
#   GITHUB_TOKEN=ghp_xxx bash push_to_your_repo.sh https://github.com/yourname/repo.git
#
# If a token is supplied (or in $GITHUB_TOKEN) it is used for this one push only
# and is NOT written to .git/config.
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_URL="${1:-}"
TOKEN="${2:-${GITHUB_TOKEN:-}}"

if [ -z "$REPO_URL" ]; then
  echo "Usage: bash push_to_your_repo.sh <repo-url> [token]" >&2
  echo "Example: bash push_to_your_repo.sh https://github.com/yourname/my-repo.git" >&2
  exit 1
fi

cd "$(dirname "$0")"
say() { printf '\n\033[1m%s\033[0m\n' "$*"; }

say "1/5  Validating results"
python3 scripts/validate_results.py | tail -3

say "2/5  Checking for copyrighted third-party text"
for bad in data/eskape_corpus.json data/eskape_triage.json data/fulltext data/eskape_fulltext; do
  [ -e "$bad" ] && { echo "ERROR: $bad present - remove before pushing" >&2; exit 1; }
done
echo "clean"

say "3/5  Initialising / updating git"
[ -d .git ] || git init -q
git symbolic-ref HEAD refs/heads/main 2>/dev/null || true
git config user.email >/dev/null 2>&1 || git config user.email "author@example.com"
git config user.name  >/dev/null 2>&1 || git config user.name  "Author"
git add -A
git diff --cached --quiet || git commit -q -m "Code and data: PAS in ESKAPE pathogens scoping review

- 419 x 31 extraction dataset (157 studies, 272 determinate observations)
- derived Table 1 values, ICC, design effect, McNemar and GEE results
- screening audit trail for 4,402 records (score + rationale, no article text)
- analysis scripts and validate_results.py (172 assertions, all passing)"

say "4/5  Configuring remote"
# strip any existing origin
git remote remove origin 2>/dev/null || true

# Build an authenticated URL only for this invocation; never persisted
PUSH_URL="$REPO_URL"
CLEAN_URL="$REPO_URL"
if [ -n "$TOKEN" ]; then
  # insert token into https URL for the push only
  PUSH_URL=$(printf '%s' "$REPO_URL" | sed -E "s#https://([^@]*@)?#https://${TOKEN}@#")
fi
git remote add origin "$CLEAN_URL"
echo "remote: $CLEAN_URL"

say "5/5  Pushing"
# push using the transient authenticated URL
git -c credential.helper= push "$PUSH_URL" HEAD:main --force

echo
say "DONE"
echo "  Repository: $CLEAN_URL"
echo
echo "Verify what you just published:"
echo "  git clone $CLEAN_URL /tmp/verify && cd /tmp/verify && python3 scripts/validate_results.py"
echo
echo "Then mint a DOI:  zenodo.org -> Settings -> GitHub -> enable repo"
echo "                  GitHub -> Releases -> tag v1.0.0"
