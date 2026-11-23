#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# publish.sh - one-shot GitHub publisher for this repository
#
# Usage:
#   bash publish.sh
#
# It will:
#   1. verify the repository passes its own validation suite
#   2. initialise git and commit
#   3. create the GitHub repository (via gh CLI if available, else prompt)
#   4. push
#
# You need: a GitHub account, and either the `gh` CLI authenticated
# (`gh auth login`) or a personal access token for the manual route.
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_NAME="${REPO_NAME:-pas-eskape-scoping-review}"
VISIBILITY="${VISIBILITY:-public}"      # public | private
BRANCH="main"

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
die() { printf '\n\033[31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
say "Step 1/5  Checking environment"
command -v python3 >/dev/null || die "python3 not found"
command -v git     >/dev/null || die "git not found"
python3 -c "import pandas, numpy, scipy" 2>/dev/null \
  || die "missing Python packages - run: pip install -r requirements.txt"

# ---------------------------------------------------------------------------
say "Step 2/5  Running the validation suite"
# The validator must pass before anything is published. A repo that does not
# reproduce its own results should never be made public.
if python3 scripts/validate_results.py > /tmp/val.log 2>&1; then
  tail -3 /tmp/val.log
else
  tail -20 /tmp/val.log
  die "validation FAILED - fix before publishing"
fi

# ---------------------------------------------------------------------------
say "Step 3/5  Checking for copyrighted third-party text"
# Guard against ever committing article abstracts/full texts.
for bad in data/eskape_corpus.json data/eskape_triage.json data/eskape_fulltext data/fulltext; do
  if [ -e "$bad" ]; then
    die "$bad exists and may contain copyrighted text. Remove it before publishing."
  fi
done
BIG=$(find . -type f -size +20M -not -path './.git/*' 2>/dev/null | head -5)
if [ -n "$BIG" ]; then
  echo "WARNING: large files found:"; echo "$BIG"
  read -r -p "Continue anyway? [y/N] " ans; [ "$ans" = "y" ] || exit 1
fi
echo "OK - no redistributable-unsafe files found"

# ---------------------------------------------------------------------------
say "Step 4/5  Initialising git"
if [ ! -d .git ]; then
  git init -q
  git branch -M "$BRANCH"
fi

# Confirm identity is set, else ask
if ! git config user.email >/dev/null 2>&1; then
  read -r -p "Your git email: " gemail;  git config user.email "$gemail"
fi
if ! git config user.name >/dev/null 2>&1; then
  read -r -p "Your name: " gname;        git config user.name  "$gname"
fi

git add .
if git diff --cached --quiet; then
  echo "nothing new to commit"
else
  git commit -q -m "PAS scoping review: extraction dataset, derived tables and analysis code

- 419 x 31 extraction dataset (157 studies, 272 determinate observations)
- derived Table 1 values, ICC, design effect, McNemar and GEE results
- screening audit trail for 4,402 records (score + rationale, no article text)
- scripts reproducing the search, screening, extraction and all statistics
- validate_results.py: 172 assertions, all passing"
  echo "committed"
fi

# ---------------------------------------------------------------------------
say "Step 5/5  Creating the GitHub repository"

if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  echo "GitHub CLI detected and authenticated."
  # Prompt for final names, defaulting to the values above
  read -r -p "Repository name [$REPO_NAME]: " rn; REPO_NAME="${rn:-$REPO_NAME}"
  read -r -p "Visibility (public/private) [$VISIBILITY]: " vis; VISIBILITY="${vis:-$VISIBILITY}"

  if git remote get-url origin >/dev/null 2>&1; then
    echo "remote 'origin' already set - pushing"
    git push -u origin "$BRANCH"
  else
    gh repo create "$REPO_NAME" --"$VISIBILITY" --source=. --remote=origin --push
  fi

  URL=$(git remote get-url origin | sed 's/\.git$//')
  say "DONE - repository published"
  echo "  URL: $URL"
  echo
  echo "Next: archive on Zenodo for a citable DOI"
  echo "  1. https://zenodo.org  ->  Settings  ->  GitHub  ->  enable this repo"
  echo "  2. GitHub  ->  Releases  ->  Create a new release  ->  tag v1.0.0"
  echo "  3. Zenodo mints a DOI automatically"
else
  echo "GitHub CLI not installed or not authenticated."
  echo
  echo "Option A (recommended) - install and authenticate gh, then re-run:"
  echo "    gh auth login"
  echo
  echo "Option B - create the repo in the browser:"
  echo "    1. Go to https://github.com/new"
  echo "    2. Repository name: $REPO_NAME"
  echo "    3. Visibility: $VISIBILITY"
  echo "    4. Do NOT add a README, .gitignore or licence (this repo has them)"
  echo "    5. Create, then run the two commands it shows you:"
  echo
  echo "         git remote add origin https://github.com/<YOUR-USERNAME>/$REPO_NAME.git"
  echo "         git push -u origin $BRANCH"
  echo
  echo "Option C - create it with a personal access token, no gh needed:"
  echo "    (create a token with 'repo' scope at https://github.com/settings/tokens)"
  echo
  echo "         export GITHUB_TOKEN=ghp_xxxxxxxxxxxx"
  echo "         curl -s -X POST https://api.github.com/user/repos \\"
  echo "              -H \"Authorization: Bearer \$GITHUB_TOKEN\" \\"
  echo "              -d '{\"name\":\"$REPO_NAME\",\"private\":false}'"
  echo "         git remote add origin https://\$GITHUB_TOKEN@github.com/<YOUR-USERNAME>/$REPO_NAME.git"
  echo "         git push -u origin $BRANCH"
fi
