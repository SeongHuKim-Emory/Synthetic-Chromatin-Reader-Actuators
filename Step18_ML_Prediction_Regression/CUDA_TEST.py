import xgboost as xgb

try:
    # Create a dummy DMatrix
    X = [[1, 2], [3, 4]]
    y = [1, 0]
    dtrain = xgb.DMatrix(X, label=y)
    
    # Attempt to train a single round on the GPU
    params = {'device': 'cuda'}
    xgb.train(params, dtrain, num_boost_round=1)
    
    print("\n✅ Success! XGBoost is correctly configured to use the GPU.")

except xgb.core.XGBoostError as e:
    if "not compiled with CUDA support" in str(e):
        print("\n❌ Failure. XGBoost is still not compiled with CUDA support.")
    else:
        print(f"\nAn unexpected XGBoost error occurred: {e}")

