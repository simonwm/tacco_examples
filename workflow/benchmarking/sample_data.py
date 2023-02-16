import tacco as tc

def sample_data(data, n_obs):
    if isinstance(n_obs, str):
        n_obs = int(n_obs)
    if n_obs > 0:
        sampled = data[tc.utils.complete_choice(data.obs.index, n_obs, seed=42)].copy()
    else:
        sampled = data.copy()
    sampled.obs_names_make_unique()
    
    return sampled