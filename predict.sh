
        TEST_FILE_PREFIX=./src/catpred_pipeline/catpred/demo/batch_kcat_input
        RECORDS_FILE=${TEST_FILE_PREFIX}.json
        CHECKPOINT_DIR=./src/catpred_pipeline/data/pretrained/production/kcat/
        
        python ./scripts/create_pdbrecords.py --data_file ${TEST_FILE_PREFIX}.csv --out_file ${RECORDS_FILE}
        python predict.py --test_path ${TEST_FILE_PREFIX}.csv --preds_path ${TEST_FILE_PREFIX}_output.csv --checkpoint_dir $CHECKPOINT_DIR --uncertainty_method mve --smiles_column SMILES --individual_ensemble_predictions --protein_records_path $RECORDS_FILE
        