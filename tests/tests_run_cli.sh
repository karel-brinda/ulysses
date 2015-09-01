#! /usr/bin/env bats

mkdir -p _output

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

@test "CLI: ulysses stats" {
	./cli/06*.sh
}

@test "CLI: ulysses query" {
	./cli/07*.sh
}

@test "CLI: ulysses shrink" {
	./cli/08*.sh
}

@test "CLI: ulysses create_many" {
	./cli/09*.sh
}
