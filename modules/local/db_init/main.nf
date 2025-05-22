process DB_INIT {
    tag "init_db"

    label 'process_low'

    input:
    val db_name

    output:
    path "${db_name}", emit: sqlite_db

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //'https://depot.galaxyproject.org/singularity/python:3.11' :
   // 'python:3.11' }"

    container 'docker.io/library/python:3.11'

    script:
    """
    mkdir -p db_init
    cat > db_init/db_setup.py << 'EOF'
import sqlite3

def init_db(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute('''CREATE TABLE IF NOT EXISTS genomes (
                      id INTEGER PRIMARY KEY,
                      genome_name TEXT NOT NULL,
                      organism TEXT)''')

    cursor.execute('''CREATE TABLE IF NOT EXISTS sketch (
                      id INTEGER PRIMARY KEY,
                      sketch_path TEXT NOT NULL)''')

    cursor.execute('''CREATE TABLE IF NOT EXISTS annotations (
                      id INTEGER PRIMARY KEY,
                      genome_id INTEGER REFERENCES genomes(id),
                      gene_name TEXT,
                      amr_annotation TEXT,
                      amr_output_path TEXT,
                      plasmid_annotation TEXT,
                      plasmid_output_path TEXT)''')

    conn.commit()
    conn.close()

if __name__ == '__main__':
    import sys
    init_db(sys.argv[1])
EOF

    python3 db_init/db_setup.py ${db_name}
    """
}
