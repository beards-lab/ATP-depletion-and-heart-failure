#!/usr/bin/env bash
set -euo pipefail

# Resolve script's own directory so it works regardless of cwd
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOKEN_FILE="$SCRIPT_DIR/data/gh_PAT.txt"
INSTALL_DIR="$HOME/.local/bin"
GIT_NAME="Filip Jezek"
GIT_EMAIL="you@example.com"   # <-- edit this once

# 1. PATH
mkdir -p "$INSTALL_DIR"
case ":$PATH:" in
  *":$INSTALL_DIR:"*) ;;
  *) export PATH="$INSTALL_DIR:$PATH"
     grep -qxF 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null \
       || echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc" ;;
esac

# 2. install gh if missing
if ! command -v gh >/dev/null 2>&1; then
  echo ">> Installing gh CLI..."
  GH_VERSION=$(curl -sIL -o /dev/null -w '%{url_effective}' \
    https://github.com/cli/cli/releases/latest \
    | grep -oE '[0-9]+\.[0-9]+\.[0-9]+$')
  echo "   latest = $GH_VERSION"
  TMP=$(mktemp -d)
  curl -fSL "https://github.com/cli/cli/releases/download/v${GH_VERSION}/gh_${GH_VERSION}_linux_amd64.tar.gz" \
    -o "$TMP/gh.tar.gz"
  tar -xzf "$TMP/gh.tar.gz" -C "$TMP"
  cp "$TMP/gh_${GH_VERSION}_linux_amd64/bin/gh" "$INSTALL_DIR/"
  chmod +x "$INSTALL_DIR/gh"
  rm -rf "$TMP"
fi
echo ">> $(gh --version | head -n1)"

# 3. auth
if [[ ! -f "$TOKEN_FILE" ]]; then
  echo "!! Token file not found: $TOKEN_FILE" >&2
  echo "   Create it: echo 'github_pat_xxx' > $TOKEN_FILE && chmod 600 $TOKEN_FILE" >&2
  exit 1
fi
chmod 600 "$TOKEN_FILE" 2>/dev/null || true
if ! gh auth status >/dev/null 2>&1; then
  echo ">> Logging in to GitHub..."
  gh auth login --with-token < "$TOKEN_FILE"
fi
gh auth status

# 4. wire into git
gh auth setup-git

# 5. identity
git config --global user.name  "$GIT_NAME"
git config --global user.email "$GIT_EMAIL"

echo ">> Done. git push should now work."
