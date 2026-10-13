# URL of the biomarker entry in the evidence source: the CIViC variant page (by CIViC variant id), or the CGI biomarker catalogue (no pages for individual CGI biomarkers)

URL of the biomarker entry in the evidence source: the CIViC variant
page (by CIViC variant id), or the CGI biomarker catalogue (no pages for
individual CGI biomarkers)

## Usage

``` r
bm_entry_url_pdf(source_db, variant_id)
```

## Arguments

- source_db:

  evidence source ('civic' or 'cgi')

- variant_id:

  variant identifier in the source (CIViC variant id)

## Value

URL per entry (NA if not available)
