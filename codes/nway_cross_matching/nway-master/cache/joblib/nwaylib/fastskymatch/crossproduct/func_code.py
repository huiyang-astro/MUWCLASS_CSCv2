# first line: 84
@mem.cache(ignore=['logger'])
def crossproduct(radectables, err, logger, pairwise_errs=[]):
	# check if away from the poles and RA=0
	use_flat_bins = True
	for ra, dec in radectables:
		if not(err < 1 and (ra > 10*err).all() and (ra < 360-10*err).all() and (numpy.abs(dec) < 45).all()):
			use_flat_bins = False
			break
	
	if use_flat_bins:
		logger.log('matching: using fast flat-sky approximation for this match')
	else:
		# choose appropriate nside for err (in deg)
		nside = 1
		for nside_next in range(30): 
			# largest distance still contained within pixels
			dist_neighbors_complete = get_healpix_resolution_degrees(2**nside_next)
			# we are looking for a pixel size which ensures bigger distances than the error radius
			# but we want the smallest pixels possible, to reduce the cartesian product
			if dist_neighbors_complete < err:
				# too small, do not accept
				# sources within err will be outside the neighbor pixels
				break
			nside = 2**nside_next
		resol = get_healpix_resolution_degrees(nside) * 60 * 60
		logger.log('matching: healpix hashing on pixel resolution ~ %f arcsec (nside=%d)' % (resol, nside))
	
	buckets = defaultdict(lambda : [[] for _ in range(len(radectables))])
	primary_cat_keys = None
	
	pbar = tqdm.tqdm(total=sum([len(t[0]) for t in radectables]))
	for ti, (ra_table, dec_table) in enumerate(radectables):
		if use_flat_bins:
			for ei, (ra, dec) in enumerate(zip(ra_table, dec_table)):
				i, j = int(ra / err), int(dec / err)
				
				# put in bucket, and neighbors
				for jj, ii in (j,i), (j,i+1), (j+1,i), (j+1, i+1):
					k = (ii, jj)
					# only primary catalogue is allowed to define new buckets
					if ti == 0 or k in buckets:
						buckets[k][ti].append(ei)
				pbar.update()
		else:
			# get healpixels
			ra, dec = ra_table, dec_table
			phi = ra / 180 * pi
			theta = dec / 180 * pi + pi/2.
			i = healpy.pixelfunc.ang2pix(nside, phi=phi, theta=theta, nest=True)
			j = healpy.pixelfunc.get_all_neighbours(nside, phi=phi, theta=theta, nest=True)
			# only consider four neighbours in one direction (N)
			# does not work, sometimes A is south of B, but B is east of A
			# so need to consider all neighbors, and deduplicate later
			neighbors = numpy.hstack((i.reshape((-1,1)), j.transpose()))
			
			# put in bucket, and neighbors
			if ti == 0:
				# only primary catalogue is allowed to define new buckets
				for ei, keys in enumerate(neighbors):
					for k in keys:
						buckets[k][ti].append(ei)
					pbar.update()
			else:
				for ei, keys in enumerate(neighbors):
					for k in keys:
						if k in primary_cat_keys:
							buckets[k][ti].append(ei)
					pbar.update()
			if ti == 0:
				primary_cat_keys = set(buckets.keys())
				
	pbar.close()
	
	# add no-counterpart options
	results = set()
	# now combine within buckets
	logger.log('matching: collecting from %d buckets, creating cartesian products ...' % len(buckets))
	#print('matching: %6d matches expected after hashing' % numpy.sum([
	#	len(lists[0]) * numpy.product([len(li) + 1 for li in lists[1:]]) 
	#		for lists in buckets.values()]))
	
	#pbar = logger.progress(ndigits=5, maxval=len(buckets)).start()
	pbar = tqdm.tqdm(total=len(buckets))
	while buckets:
		k, lists = buckets.popitem()
		pbar.update()
		# add for secondary catalogues the option of missing source
		for l in lists[1:]:
			l.insert(0, -1)
		# create the cartesian product
		local_results = itertools.product(*[sorted(l) for l in lists])
		
		# if pairwise filtering is requested, use it to trim down solutions
		if pairwise_errs:
			local_results = numpy.array(list(local_results))
			#nstart = len(local_results)
			for tablei, tablej, errij in pairwise_errs:
				indicesi = local_results[:,tablei]
				indicesj = local_results[:,tablej]
				# first find entries that actually have both entries
				mask_both = numpy.logical_and(indicesi >= 0, indicesj >= 0)
				#if not mask_both.any():
				#	continue
				# get the RA/Dec
				rai, deci = radectables[tablei]
				raj, decj = radectables[tablej]
				rai, deci = rai[indicesi[mask_both]], deci[indicesi[mask_both]]
				raj, decj = raj[indicesj[mask_both]], decj[indicesj[mask_both]]
				# compute distances
				mask_good = dist((rai, deci), (raj, decj)) < errij * 60 * 60
				# select the ones where one is missing, or those within errij
				mask_good2 = ~mask_both
				mask_good2[mask_both][mask_good] = True
				#print(mask_good2.sum(), mask_both.shape, mask_both.sum(), mask_good.sum())
				local_results = local_results[mask_good2,:]
			
			#print("compression:%d/%d" % (len(local_results), nstart))
			results.update([tuple(l) for l in local_results])
		else:
			results.update(local_results)
		del local_results
	pbar.close()

	n = len(results)
	logger.log('matching: %6d unique matches from cartesian product. sorting ...' % n)
	# now make results unique by sorting
	results = numpy.array(sorted(results))
	return results
