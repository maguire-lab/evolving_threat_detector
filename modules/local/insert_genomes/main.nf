process INSERT_GENOMES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fasta)
    path sqlite_db

    output:
    path("*.db"), emit: sqlite_db

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/python:3.11' :
        //'python:3.11' }"

    container 'docker.io/library/python:3.11'

    script:
    """
    mkdir -p scripts

    cat > scripts/insert_genomes.py << 'EOF'
import sqlite3
import sys

def insert_genome(db_path, genome_name, organism):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(
        "INSERT INTO genomes (genome_name, organism) VALUES (?, ?)",
        (genome_name, organism)
    )
    conn.commit()
    conn.close()

if __name__ == '__main__':
    db_path = sys.argv[1]
    genome_name = sys.argv[2]
    organism = sys.argv[3] if len(sys.argv) > 3 else None
    insert_genome(db_path, genome_name, organism)
EOF

    #Normalize the name to etd.db (expected by Python script)
    cp ${sqlite_db} input.db

    # Run the script
    python3 scripts/insert_genomes.py input.db ${fasta.getName()} ${meta.organism ?: ""}

    # Rename updated DB for output tracking
    cp input.db etd.db
    """
}
