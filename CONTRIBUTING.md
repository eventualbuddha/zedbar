# Contributing to zedbar

## Commit Message Format

This project uses [Conventional Commits](https://www.conventionalcommits.org/) for automated releases via release-please.

**Format**: `<type>: <description>`

### Types:
- `feat:` - New feature (triggers **minor** version bump: 0.2.0 → 0.3.0)
- `fix:` - Bug fix (triggers **patch** version bump: 0.2.0 → 0.2.1)
- `feat!:` or `fix!:` - Breaking change (triggers **major** version bump: 0.2.0 → 1.0.0)
- `chore:` - Maintenance (no release)
- `docs:` - Documentation only (no release)
- `test:` - Tests only (no release)
- `refactor:` - Code refactoring (no release)

### Examples:
```
feat: add WebP support to CLI
fix: handle empty image files correctly
feat!: remove deprecated scanImage function
chore: update dependencies
docs: add CLI usage examples
```

### Breaking Changes:
Add `!` after type OR include `BREAKING CHANGE:` in commit body:
```
feat!: redesign API interface

BREAKING CHANGE: scanImage() renamed to scanImageBytes()
```

## Release Process

1. Commit with the conventional format in a feature branch and open a PR targeting `main`
2. When PRs are merged into `main`, release-please automatically opens/updates a single release PR that includes all unreleased conventional commits since the last release
3. Merge the release PR → auto-publishes to crates.io and npm

**Note**: `feat:` and `fix:` commits are batched together into the next release based on what has been merged into `main`, rather than creating a separate release for each individual commit. Keep commits focused and atomic so the generated changelog remains clear.
