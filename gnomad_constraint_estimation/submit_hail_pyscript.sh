NAME="test"
REGION="europe-west1"
#PROJECT_ID=

hailctl dataproc start ${NAME} --region=${REGION}
gcloud config set dataproc/region ${REGION}
hailctl dataproc submit ${NAME} src/test_hailctl.py 
hailctl dataproc stop ${NAME}

