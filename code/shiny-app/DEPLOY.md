Deployment to shinyapps.io
==========================

Goal: Deploy `app.R` (this directory) to shinyapps.io including only the minimal data needed.

Why a prep step? shinyapps.io bundles ONLY the contents of the deployment directory. Your current data live outside (`../../data`). We must copy required Excel files into a local `data/` folder inside the app before deploying.

Quick Start (once per machine)

1. Install rsconnect (if not already):
   install.packages("rsconnect")

2. Get your shinyapps.io token & secret (Account > Tokens) and set them (replace placeholders):
   rsconnect::setAccountInfo(name="ACCOUNT", token="TOKEN", secret="SECRET")

Each Deployment

1. Run the prep script (copies fresh subset of data):
   source("code/shiny-app/prepare_data_for_deploy.R")

2. Deploy from the app directory (add a stable name once so future deploys update the same app):
   rsconnect::deployApp(appDir = "code/shiny-app", appName = "resazurin-explorer")

   Note: rsconnect::deployApp() does NOT have a 'title' argument. The displayed name comes from:
   - The Shiny UI (e.g. navbarPage(title = "Resazurin Trials Explorer", ...))
   - Or (optional) a DESCRIPTION file's Title field if you convert the app into a package structure.

Included Data

The prep script copies:

- trial_metadata.xlsx
- size/trial_size.xlsx
- plate-files/*/*.xlsx (EXCLUDING date 20250630 as per app logic)

Adjust what is copied by editing `prepare_data_for_deploy.R`.

Large / Sensitive Data

If data become large, instead:

- Host static .rds/.csv on GitHub (raw) or object storage (S3) and download on first run (cache locally with `tools::R_user_dir`).

Package Versions (optional but recommended)
Use renv for reproducibility:
renv::init()
renv::snapshot()
Then deploy; renv.lock will be honored by shinyapps.io.

Environment Variables (optional)
If you later externalize data, set e.g. Sys.getenv("DATA_BASE_URL") inside app.

Troubleshooting

- ERROR cannot open the connection: Confirm `data/` exists inside app directory after running prep script.
- Package not available: Make sure it's in DESCRIPTION (if you add one) or installed locally; or use renv.
- Slow startup: Pre-compute heavy summaries and save as `data/precomputed.rds`.

Contact: Add internal notes / maintainers here.
