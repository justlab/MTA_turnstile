suppressPackageStartupMessages(
   {library(data.table)
    library(stringr)
    library(pbapply)
    library(sf)})

source("../Just_universal/code/pairmemo.R")
source("../Just_universal/code/download.R")

data.root = "/data-coco/COVID_19"
pairmemo.dir = file.path(data.root, "pairmemo")
zcta.dir = "/data-belle/basemap/census/zcta"

date.first = as.Date("2014-12-27")
local.tz = "America/New_York"

crs.lonlat = 4326 # https://epsg.io/4326
crs.us.atlas = 2163 # https://epsg.io/2163

comma = scales::comma

download = function(url, to, f, ...)
    download.update.meta(data.root, url, to, f, ...)

varn = function(v, na.rm = F)
# Variance with n instead of (n - 1) in the denominator.
    mean((v - mean(v, na.rm = na.rm))^2, na.rm = na.rm)
sdn = function(v, na.rm = F)
# Standard deviation with n instead of (n - 1) in the denominator.
    sqrt(varn(v, na.rm))

turnstile = function()
  # Get the number of New York City Metropolitan Transit Authority
  # (MTA) turnstile entries and exits by station and timestamp.
  # Field descriptions:
  # https://web.archive.org/web/20190904142924/http://web.mta.info/developers/resources/nyct/turnstile/ts_Field_Description.txt
   {url.root = "http://web.mta.info/developers/data/nyct/turnstile/"
    file.interval.days = 7L

    message("Reading turnstile files")
    d = rbindlist(pblapply(
        seq(date.first, lubridate::today("America/New_York"),
            by = file.interval.days),
        function(the.date)
           {fname = sprintf("turnstile_%02d%02d%02d.txt",
              # `strftime` is platform-dependent, so it's a bit
              # safer to do this by hand.
                year(the.date) %% 1000L,
                month(the.date),
                mday(the.date))
            download(
                paste0(url.root, fname),
                file.path("mta_turnstile", fname),
                fread)}))

    message("Processing")

    setnames(d, tolower(colnames(d)))
    setnames(d, "c/a", "ca")

    stopifnot(!anyNA(d))
    stopifnot(all(d[, .(entries, exits)] >= 0))
    d = d[desc != "RECOVR AUD"]
      # Remove recovered-audit rows, which can interrupt otherwise
      # cumulative counts.
    d[, date := lubridate::mdy(date)]
    stopifnot(all(str_detect(d$time, "\\A\\d\\d:\\d\\d:\\d\\d\\z")))
    stopifnot(all(d[,
        by = .(ca, unit, scp, station, date, time), .N]$N == 1))
    d[, datetime := lubridate::ymd_hms(
        paste(date, time), tz = local.tz)]
    setkey(d, ca, unit, scp, station, datetime)

    # Observations with the same `ca` but different station names
    # should be referring to the same station. Use the latest name for
    # each station. Beware that different `ca`s can have the same
    # station name.
    message("Setting station names")
    stations = d[, by = ca,
        .(station.name = station[which.max(date)])]

    # Geolocate stations.
    locations = rbind(
        download(
            "https://raw.githubusercontent.com/chriswhong/nycturnstiles/master/geocoded.csv",
            "mta_turnstile_locations.csv",
            fread, header = F, col.names = c(
                "unit", "ca", "station", "linename", "division",
                "lat", "lon"))[!is.na(lon), .(ca, lon, lat)],
        # Add a few more station locations by hand.
        data.table(
            ca = c("N700", "N700A", "N701", "N701A", "N702", "N702A", "R550", "R551", "A049", "D001", "H007A", "N095A", "N098", "N330", "N539", "N601A", "PTH01", "PTH06", "PTH09", "PTH10", "PTH13", "PTH16", "PTH18", "PTH19", "PTH20", "PTH22", "R101", "R107D", "R108A", "R169", "R612"),
            lon = c(-73.95836, -73.95836, -73.951921, -73.951921, -73.94766, -73.94766, -74.000596, -74.000596, -74.011052, -74.011717, -73.981719, -74.007983, -74.007983, -73.86161, -73.980324, -73.966291, -74.164494, -74.033998, -74.007022, -73.998802, -73.989231, -74.027647, -74.164494, -74.164494, -74.164494, -74.011417, -74.012983, -74.011417, -74.011417, -73.972363, -73.977417),
            lat = c(40.768611, 40.768611, 40.777493, 40.777493, 40.783542, 40.783542, 40.755456, 40.755456, 40.710662, 40.635011, 40.730901, 40.709938, 40.709938, 40.729869, 40.666276, 40.764763, 40.734139, 40.726691, 40.733001, 40.734223, 40.747809, 40.735288, 40.734139, 40.734139, 40.734139, 40.71149, 40.703082, 40.71149, 40.71149, 40.79388, 40.684063)))
    stopifnot(!anyDuplicated(locations$ca))
    setkey(locations, ca)
    stations[, c("lon", "lat") :=
        locations[.(stations$ca), .(lon, lat)]]
    # Drop unlocated stations and their associated observations.
    message(sprintf("Dropping %d unlocated stations ",
        sum(is.na(stations$lon))))
    stations = stations[!is.na(lon)]
    message(sprintf("Dropping %s observations",
        d[, comma(sum(!(ca %in% stations$ca)))]))
    d = d[ca %in% stations$ca]

    # The counts of entries and exits are cumulative per (`ca`,
    # `unit`, `scp`) tuple (at least, I think that's the
    # right tuple), but sometimes the counters are reset. So, we
    # need to subtract each count from the next unless the result
    # would be negative, which indicates a reset. This method can only
    # detect resets imperfectly, but I can't see any other way to do
    # it. Worse, the number to which the counter is reset isn't
    # always 0, so we can't even assume that there have been at least
    # as many entries or exits as the new count.
    #
    # And then there are some cases that are just weird, such as
    # some turnstiles appearing to count backwards for a while,
    # (possibly due to a 64-bit integer being interpreted as 32 bits)
    # and then going back to forwards counting. We'll try to
    # exclude these cases.

    # First, remove some special cases, where counts jump back
    # and forth between observations. Probably two turnstiles
    # were being counted as one.

    d = d[!(
        (ca == "PTH12" & unit == "R542" & scp == "00-04-00" &
           date %in% (as.Date("2019-04-27") + (0:3))) |
        (ca == "PTH13" & unit == "R541" & scp == "00-00-04" &
           date %in% (as.Date("2016-02-02") + (0:2))))]

    # Now handle the rest.

    min.obs.n = 10L
    max.neg.diffs.p = .01
    min.pos.diffs.n = 3L
    max.diff = 100e3

    l = sapply(c("entries", "exits"), simplify = F, function(vname)
       {message("Decumulating - ", vname)
        n.sources = nrow(d[, by = .(ca, unit, scp), 1])
        n.sources.dropped = 0
        n.obs = nrow(d)
        n.obs.dropped = 0
        bar = txtProgressBar(style = 3, min = 0, max = n.sources)
        out = d[, by = .(ca, unit, scp),
           {setTxtProgressBar(bar, .GRP)
            x.old = get(vname)
            diffs = diff(x.old)
            if (.N >= min.obs.n &&
                    mean(diffs < 0) <= max.neg.diffs.p &&
                    sum(diffs > 0) >= min.pos.diffs.n &&
                    all(diffs <= max.diff))
               {date.diffs = as.integer(diff(date))
                i.keep = date.diffs <= 1 & diffs >= 0
                n.obs.dropped <<- n.obs.dropped + sum(!i.keep)
                .(datetime = datetime[-1][i.keep], count = diffs[i.keep])}
            else
               {n.sources.dropped <<- n.sources.dropped + 1
                n.obs.dropped <<- n.obs.dropped + .N
                NULL}}]
        close(bar)
        message(sprintf("Dropped %s sources of %s",
            comma(n.sources.dropped), comma(n.sources)))
        message(sprintf("Dropped %s observations of %s",
            comma(n.obs.dropped), comma(n.obs)))
        out[, .(ca, datetime, count)]})

    c(l, list(stations = stations))}
turnstile = pairmemo(turnstile, pairmemo.dir, mem = T)

turnstile.daily = function(hour.start = NULL, hour.end = NULL)
  # Create a data table of daily entries and exits per station.
   {l = turnstile()
    counts = sapply(c("entries", "exits"), simplify = F, function(vname)
       l[[vname]][
            (if (is.null(hour.start)) T else hour(datetime) >= hour.start) &
                (if (is.null(hour.end)) T else hour(datetime) < hour.end),
            by = .(ca, date = as.Date(datetime, tz = local.tz)),
            .(count = sum(count))])
    merge(counts[[1]], counts[[2]], all = T,
        by = c("ca", "date"),
        suffixes = paste0(".", names(counts)))}
turnstile.daily = pairmemo(turnstile.daily, pairmemo.dir, mem = T, fst = T)

station.boros = function()
  # Return a factor of boro names, with one boro for each station.
    factor(station.areas("BoroName", download(
        "https://data.cityofnewyork.us/api/geospatial/tqmj-j8zm?method=export&format=Original",
        "nyc_boros_shapefile.zip",
        function(p) read_sf(paste0("/vsizip/", p, "/nybb_20a")))))
station.boros = pairmemo(station.boros, pairmemo.dir)

station.neighborhoods = function()
  # Return an integer vector of United Hospital Fund (UHF)
  # neighborhood codes.
   {by.geo = as.integer(station.areas("UHFCODE", download(
        "https://www1.nyc.gov/assets/doh/downloads/zip/uhf42_dohmh_2009.zip",
        "nyc_uhf_nhoods_shapefile.zip",
        function(p) read_sf(paste0("/vsizip/", p, "/UHF_42_DOHMH_2009")))))
    by.hand = c(
        N037 = 304,
        N039 = 304,
        N040 = 304,
        N044 = 304,
        N045 = 304,
        N046 = 304,
        N182 = 407,
        JFK01 = 407,
        JFK02 = 407)[turnstile()$stations$ca]
    ifelse(!is.na(by.hand), as.integer(by.hand), by.geo)}
station.neighborhoods = pairmemo(station.neighborhoods, pairmemo.dir)

station.zips = function()
  # Return a vector of ZIP codes, matching locations to ZIPs on the
  # basis of ZIP Code Tabulation Areas (ZCTAs).
    {by.geo = station.areas("ZCTA5CE10", st_read(
            file.path(zcta.dir, "tl_2019_us_zcta510.shp"), quiet = T)[
        readRDS(file.path(zcta.dir, "zcta_index_in_NYC.rds")),])
     by.hand = c(
         A006 = 10022,
         A007 = 10022,
         N037 = 10025,
         N039 = 10025,
         N040 = 10025,
         N044 = 10024,
         N045 = 10024,
         N046 = 10023,
         N049 = 10019,
         N051 = 10019,
         R158 = 10019)[turnstile()$stations$ca]
     ifelse(!is.na(by.hand), as.integer(by.hand),
         as.integer(levels(by.geo))[as.integer(by.geo)])}
station.zips = pairmemo(station.zips, pairmemo.dir)

station.areas = function(colname, areas)
   {stations = st_as_sf(turnstile()$stations,
        coords = c("lon", "lat"), crs = crs.lonlat)
    result = do.call(st_intersects, lapply(list(stations, areas),
        function(x) st_transform(x, crs = crs.us.atlas)))
    areas[[colname]][sapply(result, function(i)
        if (length(i) == 0)
            NA
        else if (length(i) == 1)
            i
        else
            stop())]}

relative.subway.usage = function(the.year, by, ...)
  # For each date T in the given year, compute a measure of subway
  # usage comparing T to all other days that occur on the same month
  # and day of the week, but in different years. The usage count for T
  # is divided by the median in the other days to get a proportion.
  # The results are given per place and also pooled across places.
   {counts = copy(turnstile.daily(...))

    counts[, place := switch(by,
       boro = station.boros(),
       zcta = as.character(station.zips()),
       nhood = as.character(station.neighborhoods()),
       stop())[match(ca, turnstile()$stations$ca)]]
    # Drop observations that we can't assign a place to.
    counts = counts[!is.na(place)]

    rbindlist(lapply(c(F, T), function(place.split)
       {d = copy(counts)
        if (!place.split)
            d[, place := "all"]
        d = d[, by = .(date, place), .(uses =
            sum(count.entries, na.rm = T) +
            sum(count.exits, na.rm = T))]
        typical = d[
            year(date) != the.year,
            keyby = .(place, m = month(date), w = wday(date)),
            .(median.t = median(uses))]
        d[year(date) == the.year, .(date, place, usage =
           {x = data.table(place, month(date), wday(date))
            uses / typical[.(x), median.t]})]}))}
