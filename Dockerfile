# Used to build the test suite image. Currently, you must have the data locally, in a ./BulkData folder, for this to build successfully.
FROM rust:latest
COPY ./BulkData /Data
