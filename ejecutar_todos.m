function ejecutar_todos()
% ================================================================
%  ejecutar_todos.m
%  Ejecuta en SECUENCIA los 3 programas de la 3a aproximacion.
%  Si uno falla, captura el error, lo muestra y PASA AL SIGUIENTE.
%
%  Se estructura como FUNCION con llamadas desenrolladas (no un bucle):
%  cada programa empieza con 'clear', que borra el workspace donde corre
%  (el de 'correr'), NO el del lanzador -> la secuencia sobrevive al clear.
%
%  Uso:  >> ejecutar_todos
%  (los .m deben estar en la carpeta actual o en el path)
% ================================================================
    fprintf('\n############ LANZADOR 3a APROXIMACION ############\n');
    correr('dron_L_1D_grad_vs_ang');
    correr('dron_L_3D');
    correr('dron_escalera3D');
    fprintf('\n############ LANZADOR TERMINADO ############\n');
end

function correr(nombre)
    fprintf('\n================ EJECUTANDO: %s ================\n', nombre);
    try
        run(nombre);                       % el script corre aqui; su 'clear' borra ESTE workspace
        fprintf('\n>>> terminado SIN errores\n');   % (no se usan vars borradas por el clear)
    catch ME
        fprintf(2, '\n>>> ERROR (se continua con el siguiente):\n');
        fprintf(2, '    %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf(2, '    en %s, linea %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end
