import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// base: "./" so the built bundle works when served from any path (e.g. behind FastAPI or file://).
export default defineConfig({
  plugins: [react()],
  base: "./",
  server: { host: "127.0.0.1", port: 5178 },
  preview: { host: "127.0.0.1", port: 5178 },
});
