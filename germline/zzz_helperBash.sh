# delete the bam files which were saved with this error message:

#"{"message":"internal server error"}"
# and {"error":"You don't have access to this resource: Requested file ffaa70f9-e003-4075-af0c-c72c511d46cd does not allow read access","message":"You don't have access to this resource: Requested file ffaa70f9-e003-4075-af0c-c72c511d46cd does not allow read access"}
# they occupy less then 1Mb

find . -name "*bam" -type "f" -size -1000k -delete


