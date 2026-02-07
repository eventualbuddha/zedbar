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

1. Commit with conventional format and push to `main`
2. release-please automatically opens/updates a release PR
3. Merge the release PR → auto-publishes to crates.io and npm

**Note**: Every `feat:` or `fix:` commit that gets merged generates its own release. Keep commits focused and atomic.
