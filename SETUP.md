# NorESM / CAM sectional aerosol — repo setup on Olivia

Setup guide. It covers **two** repos:

| What | Local path | Your fork (`origin`) | Christina's (`christinafork`) |
|---|---|---|---|
| Host model + CIME | `CAM_SEC/` | `<GITHUB_USER>/CAM` | `Trudigard/CAM` |
| Sectional aerosol code (this repo) | `CAM_SEC/src/chemistry/oslo_sectional/` | `<GITHUB_USER>/OsloSectional` | `Trudigard/OsloSectional` |

## Prerequisites

- A GitHub account
- Membership of project `nn9560k` on Olivia


---

## 1. Fork on GitHub
Fork both repositories:
- **Host model:** fork  upstream (recommended)
  [`NorESMhub/CAM`](https://github.com/NorESMhub/CAM.git)  or christinas fork [`Trudigard/CAM`](https://github.com/Trudigard/CAM.git)  → `github.com/<GITHUB_USER>/CAM`
- **Aerosol code:** fork [`Trudigard/OsloSectional`](https://github.com/Trudigard/OsloSectional.git) → `github.com/<GITHUB_USER>/OsloSectional`

## 2. Clone the host model on Olivia

```bash
mkdir -p /cluster/work/projects/nn9560k/$USER
cd /cluster/work/projects/nn9560k/$USER
git clone https://github.com/<GITHUB_USER>/CAM CAM_SEC
cd CAM_SEC
```

## 3. Track Christina's fork and make a working branch

```bash
git remote add christinafork https://github.com/Trudigard/CAM.git
git fetch christinafork
git checkout -b <your_branch_name> christinafork/saltydust
```

> It's recommended to branch with `git checkout -b <name> christinafork/saltydust`. Checking out
> `christinafork/saltydust` directly leaves you in **detached HEAD**, where commits are easy to lose.

## 4. Pull in the external components
[FIRST DO THIS STEP!]( https://github.com/NorESMhub/noresm3_dev_simulations/wiki/Running-NorESM-on-Olivia)
```bash
./bin/git-fleximod update
```

This checks out `src/chemistry/oslo_sectional` and the other externals at their pinned commits.

## 5. Set up this aerosol repo for development

`git-fleximod` leaves `src/chemistry/oslo_sectional` on a detached commit pointing at Christina's
`OsloSectional`. Repoint `origin` at your fork so you can commit and push:

```bash
cd src/chemistry/oslo_sectional
git remote set-url origin https://github.com/<GITHUB_USER>/OsloSectional.git
git remote add christinafork https://github.com/Trudigard/OsloSectional.git
git fetch christinafork origin
git checkout -b wetdep_dev christinafork/saltydust
cd -
```

> After this, `git-fleximod update`
> will warn about the modified external — that is expected; do not let it reset your branch.

## 6. Run the sectional test suite



From the `CAM_SEC` root:

```bash
./cime/scripts/create_test --xml-category test_sectional --xml-machine olivia -r /cluster/work/projects/nn9560k/$USER/ -p NN9560K --output-root /cluster/work/projects/nn9560k/$USER/
```

To run a **single** test instead of the whole category, name it explicitly, e.g.:

```bash
./cime/scripts/create_test SMS_Ln9.ne16pg3_ne16pg3_mtn14.<COMPSET>.olivia_intel \
  -r /cluster/work/projects/nn9560k/$USER/ -p NN9560K \
  --output-root /cluster/work/projects/nn9560k/$USER/
```

---

## Workflow to contribute code (per now)

**Your own work**

- Commit and push to your fork / your branch — in whichever repo you changed
  (`CAM_SEC/` for host-model or CIME changes, `src/chemistry/oslo_sectional/` for aerosol code).

**Sharing upstream**

- Push your branch to your fork.
- Open a pull request against Christina's repo (`Trudigard/CAM` or `Trudigard/OsloSectional`).
