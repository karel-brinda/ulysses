#! /usr/bin/env bats

mkdir -p _output

#cd "$(dirname "$0")"
BIN=$(cd ../bin; pwd)
export PATH=$PATH:$BIN

@test "CLI: ulysses create" {
	./cli/01*.sh
}

@test "CLI: ulysses bitwise" {
	./cli/02*.sh
}

@test "CLI: ulysses hamming" {
	./cli/03*.sh
}

@test "CLI: ulysses symmdiffmat" {
	./cli/04*.sh
}

@test "CLI: ulysses dump" {
	./cli/05*.sh
}

@test "CLI: ulysses shrink" {
	./cli/06*.sh
}

@test "CLI: ulysses query" {
	./cli/07*.sh
}

@test "CLI: ulysses stats" {
	./cli/08*.sh
}

@test "CLI: ulysses create_many" {
	./cli/09*.sh
}

@test "CLI: ulysses query_and_split" {
	./cli/10*.sh
}

@test "CLI: ulysses merge" {
	./cli/11*.sh
}
