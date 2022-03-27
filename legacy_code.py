# Pasted stuff from Python

if True:
        print("Uniqueness")
        area_idx_fn = "../tmp/area_idx_500m_150vars_df.pik"
        if False:
            # prep area data
            afn = "../tmp/areaidx-osgb_grid_500m-allmat.pik"

            area_idx = _pickle_file_to_obj(afn)
            #route_idx = _pickle_file_to_obj(rfn)
            # remove columns from area
            #del route_idx
            #areas = area_idx.row_to_i.keys()
            #areas = select_random_sublist(areas, 5000)  # DEBUG
            #row_names = []
            area_df = NamedSparseMatrix.named_matrix_to_df(area_idx)
            lon_area_df = area_df[area_df.index.isin(get_london_area_ids())]
            grid_df = get_grid_geoms_df()
            lon_area_df = lon_area_df.merge(grid_df, left_index=True, right_on='id')
            lon_area_df.to_csv("../tmp/area_df_london_attr_geom.csv", index=True, index_label="area_id")
            #rows = NamedSparseMatrix.get_row_names(area_idx)
            del area_idx
            #return

            rfn = r"..\tmp\routeidx-osgb_grid_500m-all_filt_norm_df.pik"
            route_idx = _pickle_file_to_obj(rfn)
            cols_to_keep = []

            for c in route_idx.columns.values:
                if c in area_df.columns.values and not "osroutes" in c:
                    #print("\t", c)
                    cols_to_keep.append(c)
            cols_to_keep.remove("::remoteness::rem_max")
            cols_to_keep.remove("::remoteness::rem_min")
            print("Keeping columns:", len(cols_to_keep))
            area_df = area_df.loc[:, cols_to_keep]
            area_df = remove_empty_columns(area_df)
            #anorm_df = area_df
            anorm_df = transform_area_df(area_df, None)
            anorm_df.to_pickle(area_idx_fn)
            print(area_idx_fn)
            del anorm_df
            return

        if False:
            # area uniqueness
            anorm_df = _pickle_file_to_obj(area_idx_fn)
            uniq_df = calc_uniqueness_df(anorm_df, "../tmp/area_uniqueness_london", False, -1)

        if False:
            # route uniqueness
            print("route uniqueness")
            rfn = "../tmp/routeidx-osgb_grid_500m-all_filt_norm_df.pik"
            rdf = _pickle_file_to_obj(rfn)
            #rdf = rdf.sample(100, random_state=43)
            uniq_df = calc_uniqueness_df(rdf, "../tmp/route_uniqueness", True, 50000)
            return

        if True:
            print("Analyse uniqueness")
            # analyse_uniqueness("../tmp/uniqueness_area_df")
            #analyse_uniqueness("../tmp/area_uniqueness_london")
            analyse_uniqueness("../tmp/route_uniqueness")
            return
            
            
            def calc_uniqueness_df(df, fout, use_sample, sample_sz):
    print("calc_uniqueness_df", df.shape, use_sample, sample_sz, fout)
    
    class UniqParams(object):
        def uniqueness(self, row_name):
            res = {}
            sw = StopWatch(str(row_name))
            el_vec = np.reshape(self.df.loc[row_name, :].values, (-1, self.df.loc[row_name, :].size))

            if self.use_sample:
                cos_sim = sklearn.metrics.pairwise.cosine_similarity(el_vec, self.smpl_df.values)
                n = len(self.smpl_df.index)
            else:
                cos_sim = sklearn.metrics.pairwise.cosine_similarity(el_vec, self.df.values)
                n = len(self.df.index)
            uniq_idx = round(n - cos_sim.sum(), 5)
            entro = round(entropy(cos_sim[0]), 5)  # .tolist()
            #assert entro >= 0
            assert uniq_idx >= 0 and uniq_idx <= n
            res["row_name"] = row_name
            res["uniq_idx"] = round(uniq_idx / n, 5)
            res["uniq_idx_sum"] = uniq_idx
            res["vec_sum"] = round(el_vec.sum(), 5)
            res["sim_entro"] = entro
            res["sim_sum"] = round(cos_sim.sum(), 5)
            res["sim_kurto"] = round(kurtosis(cos_sim[0]), 5)
            res["sim_skew"] = round(skew(cos_sim[0]), 5)
            res["sim_min"] = round(cos_sim.min(), 5)
            res["sim_median"] = round(np.median(cos_sim[0]), 5)
            res["sim_max"] = round(cos_sim.max(), 5)
            res["is_sampled"] = self.use_sample
            # add attributes to dict
            res.update(self.df.loc[row_name, :].to_dict())
            params.counter += 1
            if params.counter % 100 == 0:
                print(params.counter)
                sw.tick()
            return res

    params = UniqParams()
    if True:
        # load row name filter
        df = df[df.index.isin(get_london_area_ids())]
        print("filtered rows", df.shape)
        if False:
            # filter columns
            col_filt = ['corine', 'restaurant', 'cafe']  # 'cafe',
            cols = df.columns
            cols_to_remove = []
            for c in cols:
                found = False
                for cc in col_filt:
                    if cc in c:
                        found = True
                if not found:
                    cols_to_remove.append(c)
            #del filt_df
            del cols, c, cc
            df = df.drop(columns=cols_to_remove)
        print(df.shape, df.columns)
        df.to_csv(fout+"_attr.csv", index=True, index_label="row_id")
    #params.uniqueness
    #params.sw = StopWatch("calc_uniqueness")
    assert len(df.index) > 0
    params.df = df
    params.use_sample = use_sample
    params.counter = 0
    if params.use_sample:
        params.smpl_df = df.sample(sample_sz, random_state=10)
        print("SAMPLE_SZ", sample_sz)

    row_names = df.index
    results = _parallel_for(row_names, 16, params.uniqueness)
    results_d = {}
    for r in results:
        results_d[r['row_name']] = r

    res_df = pd.DataFrame.from_dict(results_d, orient='index')
    res_df['uniq_rank'] = (res_df['uniq_idx']*-1).rank()
    
    fn = fout
    res_df.to_pickle(fn+'.pik')
    res_df.to_csv(fn+'.csv', index=False)
    print(res_df.shape, fn)
    return res_df


def analyse_uniqueness(fin):
    print("analyse_uniqueness", fin)
    udf = _pickle_file_to_obj(fin + ".pik")
    # density
    fn = fin + '_hist.pdf'
    udf.hist(color=None, alpha=0.5, bins=15, figsize=(30, 30))
    #plt.title("Uniqueness vars -- " + fin + " " + str(udf.shape))
    plt.savefig(fn)
    print(fn)
    plt.clf()
    plt.close()

    if "area" in fin:
        # merge area uniq with area geoms
        grid_df = get_grid_geoms_df()
        udf2 = udf.merge(grid_df, left_on='row_name', right_on='id')
        fn = fin + "_wgeoms.csv"
        print(fn, udf2.shape)
        udf2.to_csv(fn, index=False)

    if "route" in fin:
        # join route uniquness info with coords and groups
        routes_geom_df = get_route_start_pts()
        udf2 = udf.merge(routes_geom_df, left_on='row_name', right_on='route_id')

        group_df = get_route_clusters()
        print(group_df.columns, udf2.columns)
        print(group_df.shape, udf2.shape)
        print(udf.describe())
        print(udf2.describe())
        print(group_df.describe())

        udf3 = udf2.merge(group_df, left_on='row_name', right_on='route_id.1')
        fn = fin + "_groups_wgeoms.csv"
        print(fn, udf.shape, udf2.shape, udf3.shape)
        udf3.to_csv(fn, index=False)
